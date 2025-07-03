#!/usr/bin/env python3
"""
Create ATAC-seq average plots and heatmaps from BED peak files
This script creates pseudo-signal visualizations using peak data
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle
import gzip
import json

def read_bed_file(filename):
    """Read a BED file (gzipped or not) and return as DataFrame"""
    if filename.endswith('.gz'):
        with gzip.open(filename, 'rt') as f:
            df = pd.read_csv(f, sep='\t', header=None)
    else:
        df = pd.read_csv(filename, sep='\t', header=None)
    
    # Assign column names based on MACS2 BED format
    columns = ['chr', 'start', 'end', 'name', 'score', 'strand', 
               'signalValue', 'pValue', 'qValue', 'peak', 'summit', 
               'width', 'fold_enrichment']
    df.columns = columns[:len(df.columns)]
    
    return df

def create_chromosome_heatmap(g82_df, g82r_df):
    """Create a heatmap showing peak density by chromosome"""
    # Get chromosome order
    chr_order = ['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY']
    
    # Calculate peak counts and mean enrichment by chromosome
    g82_chr = g82_df.groupby('chr').agg({
        'name': 'count',
        'pValue': 'mean',
        'score': 'mean'
    }).rename(columns={'name': 'peak_count'})
    
    g82r_chr = g82r_df.groupby('chr').agg({
        'name': 'count',
        'pValue': 'mean',
        'score': 'mean'
    }).rename(columns={'name': 'peak_count'})
    
    # Create comparison matrix
    comparison_data = pd.DataFrame(index=chr_order)
    comparison_data['G82_peaks'] = g82_chr['peak_count']
    comparison_data['G82R_peaks'] = g82r_chr['peak_count']
    comparison_data['G82_enrichment'] = g82_chr['pValue']
    comparison_data['G82R_enrichment'] = g82r_chr['pValue']
    comparison_data['Peak_change'] = (comparison_data['G82R_peaks'] - comparison_data['G82_peaks']) / comparison_data['G82_peaks'] * 100
    comparison_data['Enrichment_change'] = (comparison_data['G82R_enrichment'] - comparison_data['G82_enrichment']) / comparison_data['G82_enrichment'] * 100
    
    # Filter to existing chromosomes
    comparison_data = comparison_data.dropna()
    
    # Create heatmap
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 10))
    
    # Peak count heatmap
    peak_data = comparison_data[['G82_peaks', 'G82R_peaks']].T
    sns.heatmap(peak_data, annot=True, fmt='.0f', cmap='YlOrRd', ax=ax1, cbar_kws={'label': 'Number of Peaks'})
    ax1.set_title('Peak Counts by Chromosome', fontsize=14)
    ax1.set_xlabel('Chromosome')
    ax1.set_ylabel('Sample')
    
    # Change heatmap
    change_data = comparison_data[['Peak_change', 'Enrichment_change']].T
    sns.heatmap(change_data, annot=True, fmt='.1f', cmap='RdBu_r', center=0, ax=ax2, cbar_kws={'label': 'Percentage Change (%)'})
    ax2.set_title('Percentage Changes (G82R vs G82)', fontsize=14)
    ax2.set_xlabel('Chromosome')
    ax2.set_ylabel('Metric')
    
    plt.tight_layout()
    return fig

def create_enrichment_profile(g82_df, g82r_df, top_n=1000):
    """Create average enrichment profile for top peaks"""
    # Get top peaks by enrichment score
    g82_top = g82_df.nlargest(top_n, 'pValue').sort_values(['chr', 'start'])
    g82r_top = g82r_df.nlargest(top_n, 'pValue').sort_values(['chr', 'start'])
    
    # Create profile plot
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
    
    # Plot 1: Average enrichment by peak rank
    window = 50  # Moving average window
    g82_sorted = g82_df.sort_values('pValue', ascending=False)['pValue'].values
    g82r_sorted = g82r_df.sort_values('pValue', ascending=False)['pValue'].values
    
    # Calculate moving averages
    n_peaks = min(10000, len(g82_sorted), len(g82r_sorted))
    x = np.arange(n_peaks)
    
    g82_ma = pd.Series(g82_sorted[:n_peaks]).rolling(window=window, center=True).mean()
    g82r_ma = pd.Series(g82r_sorted[:n_peaks]).rolling(window=window, center=True).mean()
    
    ax1.plot(x, g82_ma, label='G82 (Parental)', color='#3498db', linewidth=2)
    ax1.plot(x, g82r_ma, label='G82R (Resistant)', color='#e74c3c', linewidth=2)
    ax1.set_xlabel('Peak Rank')
    ax1.set_ylabel('Enrichment Score (-log10 p-value)')
    ax1.set_title(f'Average Enrichment Profile (Top {n_peaks} peaks)', fontsize=14)
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Peak width vs enrichment
    ax2.scatter(g82_df['width'], g82_df['pValue'], alpha=0.3, s=5, label='G82', color='#3498db')
    ax2.scatter(g82r_df['width'], g82r_df['pValue'], alpha=0.3, s=5, label='G82R', color='#e74c3c')
    ax2.set_xlabel('Peak Width (bp)')
    ax2.set_ylabel('Enrichment Score (-log10 p-value)')
    ax2.set_title('Peak Width vs Enrichment Score', fontsize=14)
    ax2.legend()
    ax2.set_xlim(0, 2000)  # Focus on typical peak widths
    ax2.set_ylim(0, 50)    # Focus on reasonable enrichment range
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    return fig

def create_signal_heatmap(g82_df, g82r_df, n_regions=500):
    """Create a heatmap of signal intensities for top differential regions"""
    # Find differential regions (peaks with largest fold change in signal)
    g82_peaks = g82_df.set_index(['chr', 'start', 'end'])
    g82r_peaks = g82r_df.set_index(['chr', 'start', 'end'])
    
    # Get common peaks (approximate overlap)
    common_peaks = []
    for idx in g82_peaks.index:
        chr_match = g82r_peaks[g82r_peaks.index.get_level_values('chr') == idx[0]]
        if not chr_match.empty:
            # Find overlapping peaks
            for r_idx in chr_match.index:
                if (r_idx[1] <= idx[2]) and (r_idx[2] >= idx[1]):  # Overlap condition
                    common_peaks.append({
                        'chr': idx[0],
                        'start': max(idx[1], r_idx[1]),
                        'end': min(idx[2], r_idx[2]),
                        'g82_signal': g82_peaks.loc[idx, 'score'],
                        'g82r_signal': g82r_peaks.loc[r_idx, 'score'],
                        'g82_enrichment': g82_peaks.loc[idx, 'pValue'],
                        'g82r_enrichment': g82r_peaks.loc[r_idx, 'pValue']
                    })
    
    common_df = pd.DataFrame(common_peaks)
    if not common_df.empty:
        common_df['signal_fc'] = np.log2((common_df['g82r_signal'] + 1) / (common_df['g82_signal'] + 1))
        common_df = common_df.sort_values('signal_fc', key=abs, ascending=False).head(n_regions)
    
    # Create heatmap data
    fig, ax = plt.subplots(figsize=(10, 12))
    
    if not common_df.empty:
        # Prepare data for heatmap
        heatmap_data = common_df[['g82_signal', 'g82r_signal']].T
        heatmap_data.index = ['G82 (Parental)', 'G82R (Resistant)']
        
        # Create heatmap
        sns.heatmap(heatmap_data, cmap='YlOrRd', cbar_kws={'label': 'Signal Intensity'}, ax=ax)
        ax.set_title(f'Signal Intensity Heatmap (Top {len(common_df)} Differential Regions)', fontsize=14)
        ax.set_xlabel('Genomic Regions (sorted by fold change)')
        ax.set_ylabel('Sample')
        
        # Don't show all x-tick labels if too many
        if len(common_df) > 50:
            ax.set_xticks([])
    else:
        ax.text(0.5, 0.5, 'No overlapping peaks found for comparison', 
                ha='center', va='center', transform=ax.transAxes, fontsize=14)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
    
    plt.tight_layout()
    return fig

def create_genomic_distribution_plot(g82_df, g82r_df):
    """Create a plot showing distribution of peaks across genome"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Plot 1: Peak density along chromosomes
    chr_order = ['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY']
    chr_lengths = {  # Approximate human chromosome lengths (hg38)
        'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559, 'chr4': 190214555,
        'chr5': 181538259, 'chr6': 170805979, 'chr7': 159345973, 'chr8': 145138636,
        'chr9': 138394717, 'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309,
        'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189, 'chr16': 90338345,
        'chr17': 83257441, 'chr18': 80373285, 'chr19': 58617616, 'chr20': 64444167,
        'chr21': 46709983, 'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415
    }
    
    # Calculate peak density (peaks per Mb)
    g82_density = g82_df.groupby('chr').size() / pd.Series(chr_lengths) * 1e6
    g82r_density = g82r_df.groupby('chr').size() / pd.Series(chr_lengths) * 1e6
    
    x = np.arange(len(chr_order))
    width = 0.35
    
    g82_vals = [g82_density.get(chr, 0) for chr in chr_order]
    g82r_vals = [g82r_density.get(chr, 0) for chr in chr_order]
    
    ax1.bar(x - width/2, g82_vals, width, label='G82', color='#3498db')
    ax1.bar(x + width/2, g82r_vals, width, label='G82R', color='#e74c3c')
    ax1.set_xlabel('Chromosome')
    ax1.set_ylabel('Peak Density (peaks per Mb)')
    ax1.set_title('Peak Density by Chromosome', fontsize=14)
    ax1.set_xticks(x)
    ax1.set_xticklabels(chr_order, rotation=45, ha='right')
    ax1.legend()
    ax1.grid(True, alpha=0.3, axis='y')
    
    # Plot 2: Peak size distribution
    ax2.hist(g82_df['width'], bins=50, alpha=0.7, label='G82', color='#3498db', density=True)
    ax2.hist(g82r_df['width'], bins=50, alpha=0.7, label='G82R', color='#e74c3c', density=True)
    ax2.set_xlabel('Peak Width (bp)')
    ax2.set_ylabel('Density')
    ax2.set_title('Peak Width Distribution', fontsize=14)
    ax2.set_xlim(0, 2000)
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    return fig

def main():
    # Read BED files
    print("Reading BED files...")
    g82_df = read_bed_file('GSM3746025_G82.bed.gz')
    g82r_df = read_bed_file('GSM3746026_G82R.bed.gz')
    
    print(f"G82 peaks: {len(g82_df)}")
    print(f"G82R peaks: {len(g82r_df)}")
    
    # Create visualizations
    print("Creating chromosome heatmap...")
    fig1 = create_chromosome_heatmap(g82_df, g82r_df)
    fig1.savefig('atac_chromosome_heatmap.png', dpi=300, bbox_inches='tight')
    
    print("Creating enrichment profile...")
    fig2 = create_enrichment_profile(g82_df, g82r_df)
    fig2.savefig('atac_enrichment_profile.png', dpi=300, bbox_inches='tight')
    
    print("Creating signal heatmap...")
    fig3 = create_signal_heatmap(g82_df, g82r_df)
    fig3.savefig('atac_signal_heatmap.png', dpi=300, bbox_inches='tight')
    
    print("Creating genomic distribution plot...")
    fig4 = create_genomic_distribution_plot(g82_df, g82r_df)
    fig4.savefig('atac_genomic_distribution.png', dpi=300, bbox_inches='tight')
    
    print("All visualizations created successfully!")
    
    # Save plot data for HTML integration
    plot_data = {
        'chromosome_heatmap': 'atac_chromosome_heatmap.png',
        'enrichment_profile': 'atac_enrichment_profile.png',
        'signal_heatmap': 'atac_signal_heatmap.png',
        'genomic_distribution': 'atac_genomic_distribution.png'
    }
    
    with open('atac_plot_paths.json', 'w') as f:
        json.dump(plot_data, f, indent=2)

if __name__ == '__main__':
    main()