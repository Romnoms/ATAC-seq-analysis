#!/usr/bin/env python3
"""
ATAC-seq Peak Visualization Script
Generates comprehensive plots for ATAC-seq peak analysis
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib_venn import venn2
import gzip
import argparse
import os
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Set publication-quality style
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10
plt.rcParams['font.family'] = 'sans-serif'


def read_bed_file(filename):
    """Read MACS2 narrowPeak BED file format"""
    with gzip.open(filename, 'rt') as f:
        data = []
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 10:
                data.append({
                    'chr': parts[0],
                    'start': int(parts[1]),
                    'end': int(parts[2]),
                    'name': parts[3],
                    'score': int(parts[4]),
                    'strand': parts[5] if len(parts) > 5 else '.',
                    'signalValue': float(parts[6]) if len(parts) > 6 else 0,
                    'pValue': float(parts[7]) if len(parts) > 7 else 0,
                    'qValue': float(parts[8]) if len(parts) > 8 else 0,
                    'peak': int(parts[9]) if len(parts) > 9 else 0
                })
    
    df = pd.DataFrame(data)
    df['width'] = df['end'] - df['start']
    return df


def plot_peak_counts(g82_data, g82r_data, output_dir):
    """Generate peak count comparison bar plot"""
    fig, ax = plt.subplots(figsize=(8, 6))
    
    samples = ['G82\n(Parental)', 'G82R\n(Resistant)']
    counts = [len(g82_data), len(g82r_data)]
    colors = ['#3498db', '#e74c3c']
    
    bars = ax.bar(samples, counts, color=colors, width=0.6, edgecolor='black', linewidth=1.5)
    
    # Add value labels
    for bar, count in zip(bars, counts):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 500,
                f'{count:,}', ha='center', va='bottom', fontsize=12, fontweight='bold')
    
    # Calculate percentage change
    pct_change = ((counts[1] - counts[0]) / counts[0]) * 100
    ax.text(0.5, 0.95, f'Change: +{pct_change:.1f}%', transform=ax.transAxes,
            ha='center', va='top', fontsize=11, style='italic',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='yellow', alpha=0.3))
    
    ax.set_ylabel('Number of Peaks', fontsize=14)
    ax.set_title('Total ATAC-seq Peaks', fontsize=16, fontweight='bold')
    ax.set_ylim(0, max(counts) * 1.15)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'peak_count_comparison.png'), bbox_inches='tight')
    plt.close()


def plot_peak_width_distribution(g82_data, g82r_data, output_dir):
    """Generate peak width distribution histogram"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Histogram
    ax1.hist(g82_data['width'], bins=50, alpha=0.7, label='G82 (Parental)', 
             color='#3498db', range=(0, 1500), density=True)
    ax1.hist(g82r_data['width'], bins=50, alpha=0.7, label='G82R (Resistant)', 
             color='#e74c3c', range=(0, 1500), density=True)
    
    ax1.set_xlabel('Peak Width (bp)', fontsize=12)
    ax1.set_ylabel('Density', fontsize=12)
    ax1.set_title('Peak Width Distribution', fontsize=14, fontweight='bold')
    ax1.legend()
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    
    # Box plot
    box_data = pd.DataFrame({
        'Width': pd.concat([g82_data['width'], g82r_data['width']]),
        'Sample': ['G82'] * len(g82_data) + ['G82R'] * len(g82r_data)
    })
    
    sns.boxplot(x='Sample', y='Width', data=box_data, ax=ax2, 
                palette=['#3498db', '#e74c3c'])
    ax2.set_ylabel('Peak Width (bp)', fontsize=12)
    ax2.set_title('Peak Width Comparison', fontsize=14, fontweight='bold')
    ax2.set_ylim(0, 1500)
    
    # Add mean values
    means = box_data.groupby('Sample')['Width'].mean()
    for i, (sample, mean) in enumerate(means.items()):
        ax2.text(i, 1400, f'Mean: {mean:.0f}', ha='center', fontsize=10, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'peak_width_distribution.png'), bbox_inches='tight')
    plt.close()


def plot_chromosome_changes(g82_data, g82r_data, output_dir):
    """Generate chromosome-specific peak change plot"""
    # Standard chromosome order
    chr_order = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
    
    # Count peaks per chromosome
    g82_chr = g82_data['chr'].value_counts()
    g82r_chr = g82r_data['chr'].value_counts()
    
    # Calculate changes
    chr_changes = []
    for chr in chr_order:
        if chr in g82_chr.index or chr in g82r_chr.index:
            g82_count = g82_chr.get(chr, 0)
            g82r_count = g82r_chr.get(chr, 0)
            if g82_count > 0:
                pct_change = ((g82r_count - g82_count) / g82_count) * 100
            else:
                pct_change = 100 if g82r_count > 0 else 0
            
            chr_changes.append({
                'chr': chr,
                'pct_change': pct_change,
                'g82': g82_count,
                'g82r': g82r_count,
                'diff': g82r_count - g82_count
            })
    
    chr_df = pd.DataFrame(chr_changes)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Color based on change direction
    colors = ['#e74c3c' if x < 0 else '#27ae60' for x in chr_df['pct_change']]
    
    bars = ax.barh(chr_df['chr'], chr_df['pct_change'], color=colors, edgecolor='black', linewidth=0.5)
    ax.axvline(x=0, color='black', linestyle='--', alpha=0.5)
    
    # Add value labels for significant changes
    for i, row in chr_df.iterrows():
        if abs(row['pct_change']) > 10:
            x_pos = row['pct_change'] + (2 if row['pct_change'] > 0 else -2)
            ax.text(x_pos, i, f"{row['diff']:+d}", ha='left' if row['pct_change'] > 0 else 'right',
                   va='center', fontsize=9, fontweight='bold')
    
    ax.set_xlabel('Percentage Change (%)', fontsize=14)
    ax.set_ylabel('Chromosome', fontsize=14)
    ax.set_title('Chromosome-specific Peak Changes', fontsize=16, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlim(-30, 30)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'chromosome_peak_changes.png'), bbox_inches='tight')
    plt.close()


def plot_signal_distribution(g82_data, g82r_data, output_dir):
    """Generate signal intensity distribution plots"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Filter extreme outliers for visualization
    g82_signal = g82_data['score'][g82_data['score'] < g82_data['score'].quantile(0.99)]
    g82r_signal = g82r_data['score'][g82r_data['score'] < g82r_data['score'].quantile(0.99)]
    
    # Density plot
    ax1.hist(g82_signal, bins=50, alpha=0.7, density=True, label='G82 (Parental)', color='#3498db')
    ax1.hist(g82r_signal, bins=50, alpha=0.7, density=True, label='G82R (Resistant)', color='#e74c3c')
    
    ax1.set_xlabel('Signal Intensity (Score)', fontsize=12)
    ax1.set_ylabel('Density', fontsize=12)
    ax1.set_title('Signal Intensity Distribution', fontsize=14, fontweight='bold')
    ax1.legend()
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    
    # Violin plot for -log10(p-value)
    pval_data = pd.DataFrame({
        'pValue': pd.concat([g82_data['pValue'], g82r_data['pValue']]),
        'Sample': ['G82'] * len(g82_data) + ['G82R'] * len(g82r_data)
    })
    
    # Cap extreme values for visualization
    pval_data['pValue'] = pval_data['pValue'].clip(upper=100)
    
    sns.violinplot(x='Sample', y='pValue', data=pval_data, ax=ax2,
                   palette=['#3498db', '#e74c3c'])
    ax2.set_ylabel('Enrichment Score (-log10 p-value)', fontsize=12)
    ax2.set_title('Peak Enrichment Distribution', fontsize=14, fontweight='bold')
    ax2.set_ylim(0, 100)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'signal_distribution.png'), bbox_inches='tight')
    plt.close()


def plot_peak_overlap(g82_data, g82r_data, output_dir):
    """Generate Venn diagram for peak overlap"""
    from matplotlib_venn import venn2
    
    # Create peak identifiers (using a window for overlap)
    def get_peak_regions(df, window=100):
        regions = set()
        for _, row in df.iterrows():
            # Round to nearest window size for approximate overlap
            start = (row['start'] // window) * window
            end = ((row['end'] + window - 1) // window) * window
            regions.add(f"{row['chr']}:{start}-{end}")
        return regions
    
    g82_regions = get_peak_regions(g82_data)
    g82r_regions = get_peak_regions(g82r_data)
    
    fig, ax = plt.subplots(figsize=(8, 8))
    
    venn = venn2([g82_regions, g82r_regions], 
                 set_labels=('G82 (Parental)', 'G82R (Resistant)'),
                 set_colors=('#3498db', '#e74c3c'), 
                 alpha=0.7)
    
    # Customize labels
    for text in venn.set_labels:
        text.set_fontsize(14)
        text.set_fontweight('bold')
    
    for text in venn.subset_labels:
        if text:
            text.set_fontsize(12)
    
    # Calculate overlap percentage
    overlap = len(g82_regions.intersection(g82r_regions))
    total = len(g82_regions.union(g82r_regions))
    overlap_pct = (overlap / len(g82r_regions)) * 100
    
    ax.set_title('Peak Overlap Analysis', fontsize=16, fontweight='bold')
    ax.text(0.5, 0.02, f'Overlap: {overlap_pct:.1f}% of G82R peaks', 
            transform=ax.transAxes, ha='center', fontsize=11, style='italic')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'peak_overlap_venn.png'), bbox_inches='tight')
    plt.close()


def generate_summary_table(g82_data, g82r_data, output_dir):
    """Generate summary statistics table"""
    summary = {
        'Metric': ['Total Peaks', 'Mean Width (bp)', 'Median Width (bp)', 
                   'Mean Score', 'Mean -log10(p-value)', 'Peaks > 500bp', 
                   'Peaks > 1000bp'],
        'G82 (Parental)': [
            len(g82_data),
            f"{g82_data['width'].mean():.1f}",
            f"{g82_data['width'].median():.0f}",
            f"{g82_data['score'].mean():.1f}",
            f"{g82_data['pValue'].mean():.1f}",
            len(g82_data[g82_data['width'] > 500]),
            len(g82_data[g82_data['width'] > 1000])
        ],
        'G82R (Resistant)': [
            len(g82r_data),
            f"{g82r_data['width'].mean():.1f}",
            f"{g82r_data['width'].median():.0f}",
            f"{g82r_data['score'].mean():.1f}",
            f"{g82r_data['pValue'].mean():.1f}",
            len(g82r_data[g82r_data['width'] > 500]),
            len(g82r_data[g82r_data['width'] > 1000])
        ]
    }
    
    summary_df = pd.DataFrame(summary)
    
    # Calculate percentage changes
    changes = []
    for i in range(len(summary_df)):
        if i == 0:  # Total peaks
            val1 = len(g82_data)
            val2 = len(g82r_data)
        elif i in [1, 2, 3, 4]:  # Means
            val1 = float(summary_df.iloc[i, 1].replace(',', ''))
            val2 = float(summary_df.iloc[i, 2].replace(',', ''))
        else:  # Counts
            val1 = int(summary_df.iloc[i, 1])
            val2 = int(summary_df.iloc[i, 2])
        
        if val1 > 0:
            pct_change = ((val2 - val1) / val1) * 100
            changes.append(f"{pct_change:+.1f}%")
        else:
            changes.append("N/A")
    
    summary_df['Change (%)'] = changes
    
    # Save as CSV
    summary_df.to_csv(os.path.join(output_dir, 'summary_statistics.csv'), index=False)
    
    return summary_df


def main():
    parser = argparse.ArgumentParser(description='Generate ATAC-seq peak analysis plots')
    parser.add_argument('--g82', type=str, 
                       default='data/GSM3746025_G82.bed.gz',
                       help='Path to G82 (parental) BED file')
    parser.add_argument('--g82r', type=str,
                       default='data/GSM3746026_G82R.bed.gz',
                       help='Path to G82R (resistant) BED file')
    parser.add_argument('--output', type=str, default='results/figures',
                       help='Output directory for plots')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output, exist_ok=True)
    
    print("ATAC-seq Peak Analysis")
    print("=" * 50)
    
    # Read data
    print(f"Reading G82 data from {args.g82}...")
    g82_data = read_bed_file(args.g82)
    print(f"  - Loaded {len(g82_data):,} peaks")
    
    print(f"Reading G82R data from {args.g82r}...")
    g82r_data = read_bed_file(args.g82r)
    print(f"  - Loaded {len(g82r_data):,} peaks")
    
    # Generate plots
    print("\nGenerating visualizations...")
    
    print("  - Peak count comparison")
    plot_peak_counts(g82_data, g82r_data, args.output)
    
    print("  - Peak width distribution")
    plot_peak_width_distribution(g82_data, g82r_data, args.output)
    
    print("  - Chromosome-specific changes")
    plot_chromosome_changes(g82_data, g82r_data, args.output)
    
    print("  - Signal intensity distribution")
    plot_signal_distribution(g82_data, g82r_data, args.output)
    
    print("  - Peak overlap analysis")
    plot_peak_overlap(g82_data, g82r_data, args.output)
    
    print("  - Summary statistics")
    summary_df = generate_summary_table(g82_data, g82r_data, args.output)
    
    print("\nSummary Statistics:")
    print(summary_df.to_string(index=False))
    
    print(f"\nAll plots saved to: {args.output}/")
    print("Analysis complete!")


if __name__ == "__main__":
    main()