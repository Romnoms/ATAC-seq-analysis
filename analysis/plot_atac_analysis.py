#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib_venn import venn2
import gzip
import warnings
warnings.filterwarnings('ignore')

# Set style
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")

# Create output directory
import os
os.makedirs('plots', exist_ok=True)

# Read the BED files
def read_bed_file(filename):
    with gzip.open(filename, 'rt') as f:
        data = []
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 13:
                data.append({
                    'chr': parts[0],
                    'start': int(parts[1]),
                    'end': int(parts[2]),
                    'name': parts[3],
                    'score': int(parts[4]),
                    'strand': parts[5],
                    'signalValue': float(parts[6]),
                    'pValue': float(parts[7]),
                    'qValue': float(parts[8]),
                    'peak': int(parts[9]),
                    'summit': int(parts[10]),
                    'width': int(parts[11]),
                    'fold_enrichment': float(parts[12])
                })
    return pd.DataFrame(data)

print("Reading BED files...")
g82_data = read_bed_file('GSM3746025_G82.bed.gz')
g82r_data = read_bed_file('GSM3746026_G82R.bed.gz')

# Calculate peak widths
g82_data['peak_width'] = g82_data['end'] - g82_data['start']
g82r_data['peak_width'] = g82r_data['end'] - g82r_data['start']

# 1. Peak Count Comparison Bar Plot
fig, ax = plt.subplots(figsize=(8, 6))
samples = ['G82\n(Parental)', 'G82R\n(Resistant)']
counts = [len(g82_data), len(g82r_data)]
colors = ['#3498db', '#e74c3c']

bars = ax.bar(samples, counts, color=colors, width=0.6)
for bar, count in zip(bars, counts):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 500,
            f'{count:,}', ha='center', va='bottom', fontsize=12, fontweight='bold')

ax.set_ylabel('Number of Peaks', fontsize=14)
ax.set_title('Total ATAC-seq Peaks', fontsize=16, fontweight='bold')
ax.set_subtitle('Comparison between Parental and Resistant Cells', fontsize=12)
ax.set_ylim(0, max(counts) * 1.1)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig('plots/01_peak_count_comparison.png', dpi=300, bbox_inches='tight')
plt.close()

# 2. Peak Width Distribution
fig, ax = plt.subplots(figsize=(10, 6))
ax.hist(g82_data['peak_width'], bins=50, alpha=0.7, label='G82 (Parental)', 
        color='#3498db', range=(0, 1500))
ax.hist(g82r_data['peak_width'], bins=50, alpha=0.7, label='G82R (Resistant)', 
        color='#e74c3c', range=(0, 1500))

ax.set_xlabel('Peak Width (bp)', fontsize=14)
ax.set_ylabel('Frequency', fontsize=14)
ax.set_title('Peak Width Distribution', fontsize=16, fontweight='bold')
ax.set_subtitle('Comparison of ATAC-seq Peak Widths', fontsize=12)
ax.legend()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig('plots/02_peak_width_distribution.png', dpi=300, bbox_inches='tight')
plt.close()

# 3. Chromosome Distribution
# Count peaks per chromosome
chr_order = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
g82_chr_counts = g82_data['chr'].value_counts()
g82r_chr_counts = g82r_data['chr'].value_counts()

# Calculate percentage change
chr_changes = []
for chr in chr_order:
    if chr in g82_chr_counts and chr in g82r_chr_counts:
        g82_count = g82_chr_counts[chr]
        g82r_count = g82r_chr_counts[chr]
        pct_change = ((g82r_count - g82_count) / g82_count) * 100
        chr_changes.append({'chr': chr, 'pct_change': pct_change, 
                           'g82': g82_count, 'g82r': g82r_count})

chr_df = pd.DataFrame(chr_changes)

fig, ax = plt.subplots(figsize=(10, 8))
colors = ['#e74c3c' if x < 0 else '#27ae60' for x in chr_df['pct_change']]
bars = ax.barh(chr_df['chr'], chr_df['pct_change'], color=colors)
ax.axvline(x=0, color='black', linestyle='--', alpha=0.5)

ax.set_xlabel('Percentage Change (%)', fontsize=14)
ax.set_ylabel('Chromosome', fontsize=14)
ax.set_title('Chromosome-specific Peak Changes', fontsize=16, fontweight='bold')
ax.set_subtitle('Percentage change in peak counts (G82R vs G82)', fontsize=12)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig('plots/03_chromosome_peak_changes.png', dpi=300, bbox_inches='tight')
plt.close()

# 4. Enrichment Score Distribution
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

# Box plot
box_data = [g82_data['pValue'], g82r_data['pValue']]
bp = ax1.boxplot(box_data, labels=['G82', 'G82R'], patch_artist=True, showfliers=False)
bp['boxes'][0].set_facecolor('#3498db')
bp['boxes'][1].set_facecolor('#e74c3c')

ax1.set_ylabel('Enrichment Score (-log10 p-value)', fontsize=12)
ax1.set_title('Box Plot', fontsize=14)
ax1.set_ylim(0, 100)

# Violin plot
parts = ax2.violinplot(box_data, positions=[1, 2], showmeans=True, showextrema=True)
parts['bodies'][0].set_facecolor('#3498db')
parts['bodies'][1].set_facecolor('#e74c3c')

ax2.set_xticks([1, 2])
ax2.set_xticklabels(['G82', 'G82R'])
ax2.set_ylabel('Enrichment Score (-log10 p-value)', fontsize=12)
ax2.set_title('Violin Plot', fontsize=14)
ax2.set_ylim(0, 100)

fig.suptitle('Enrichment Score Distribution', fontsize=16, fontweight='bold')
plt.tight_layout()
plt.savefig('plots/04_enrichment_distribution.png', dpi=300, bbox_inches='tight')
plt.close()

# 5. Signal Intensity Distribution
fig, ax = plt.subplots(figsize=(10, 6))

# Filter outliers for better visualization
g82_signal = g82_data['score'][g82_data['score'] < 1000]
g82r_signal = g82r_data['score'][g82r_data['score'] < 1000]

ax.hist(g82_signal, bins=50, alpha=0.7, density=True, label='G82 (Parental)', color='#3498db')
ax.hist(g82r_signal, bins=50, alpha=0.7, density=True, label='G82R (Resistant)', color='#e74c3c')

ax.set_xlabel('Signal Intensity', fontsize=14)
ax.set_ylabel('Density', fontsize=14)
ax.set_title('Signal Intensity Distribution', fontsize=16, fontweight='bold')
ax.set_subtitle('ATAC-seq Peak Signal Scores', fontsize=12)
ax.legend()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig('plots/05_signal_intensity_distribution.png', dpi=300, bbox_inches='tight')
plt.close()

# 6. Venn Diagram for Peak Overlap
# Create peak identifiers
g82_peaks = set(f"{row['chr']}:{row['start']}-{row['end']}" for _, row in g82_data.iterrows())
g82r_peaks = set(f"{row['chr']}:{row['start']}-{row['end']}" for _, row in g82r_data.iterrows())

# Calculate overlap (approximate)
overlap = len(g82_peaks.intersection(g82r_peaks))

fig, ax = plt.subplots(figsize=(8, 8))
venn = venn2([g82_peaks, g82r_peaks], set_labels=('G82 (Parental)', 'G82R (Resistant)'),
             set_colors=('#3498db', '#e74c3c'), alpha=0.7)

# Customize labels
for text in venn.set_labels:
    text.set_fontsize(14)
for text in venn.subset_labels:
    if text:
        text.set_fontsize(12)

ax.set_title('Peak Overlap Analysis', fontsize=16, fontweight='bold')
ax.set_subtitle('Shared and unique ATAC-seq peaks', fontsize=12)

plt.tight_layout()
plt.savefig('plots/06_peak_overlap_venn.png', dpi=300, bbox_inches='tight')
plt.close()

# 7. Summary Statistics Heatmap
summary_data = {
    'Total Peaks': [len(g82_data), len(g82r_data)],
    'Mean Width': [g82_data['peak_width'].mean(), g82r_data['peak_width'].mean()],
    'Mean Enrichment': [g82_data['pValue'].mean(), g82r_data['pValue'].mean()],
    'Mean Signal': [g82_data['score'].mean(), g82r_data['score'].mean()]
}

summary_df = pd.DataFrame(summary_data, index=['G82', 'G82R'])

# Normalize for heatmap
summary_norm = summary_df.copy()
for col in summary_norm.columns:
    min_val = summary_norm[col].min()
    max_val = summary_norm[col].max()
    summary_norm[col] = (summary_norm[col] - min_val) / (max_val - min_val)

fig, ax = plt.subplots(figsize=(10, 4))
im = ax.imshow(summary_norm.values, cmap='RdBu_r', aspect='auto')

# Set ticks and labels
ax.set_xticks(np.arange(len(summary_df.columns)))
ax.set_yticks(np.arange(len(summary_df.index)))
ax.set_xticklabels(summary_df.columns)
ax.set_yticklabels(summary_df.index)

# Add text annotations
for i in range(len(summary_df.index)):
    for j in range(len(summary_df.columns)):
        text = ax.text(j, i, f'{summary_df.iloc[i, j]:.1f}',
                      ha='center', va='center', color='white', fontweight='bold')

ax.set_title('Summary Statistics Comparison', fontsize=16, fontweight='bold')
ax.set_subtitle('Key metrics between G82 and G82R samples', fontsize=12)

# Add colorbar
cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
cbar.set_label('Relative Value', fontsize=12)

plt.tight_layout()
plt.savefig('plots/07_summary_statistics_heatmap.png', dpi=300, bbox_inches='tight')
plt.close()

# 8. Top Peaks Comparison
top_n = 20
g82_top = g82_data.nlargest(top_n, 'pValue').reset_index(drop=True)
g82r_top = g82r_data.nlargest(top_n, 'pValue').reset_index(drop=True)

fig, ax = plt.subplots(figsize=(10, 6))
x = np.arange(1, top_n + 1)
ax.plot(x, g82_top['pValue'], 'o-', color='#3498db', linewidth=2, 
        markersize=8, label='G82 (Parental)')
ax.plot(x, g82r_top['pValue'], 's-', color='#e74c3c', linewidth=2, 
        markersize=8, label='G82R (Resistant)')

ax.set_xlabel('Peak Rank', fontsize=14)
ax.set_ylabel('Enrichment Score (-log10 p-value)', fontsize=14)
ax.set_title('Top 20 Peaks by Enrichment Score', fontsize=16, fontweight='bold')
ax.set_subtitle('Comparison of highest scoring peaks', fontsize=12)
ax.legend()
ax.grid(True, alpha=0.3)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig('plots/08_top_peaks_comparison.png', dpi=300, bbox_inches='tight')
plt.close()

# Print summary
print("\nATAC-seq Analysis Plots Generated Successfully!")
print("===============================================")
print("Output directory: ./plots/")
print("\nGenerated plots:")
print("1. Peak count comparison")
print("2. Peak width distribution")
print("3. Chromosome-specific changes")
print("4. Enrichment score distribution")
print("5. Signal intensity distribution")
print("6. Peak overlap Venn diagram")
print("7. Summary statistics heatmap")
print("8. Top peaks comparison")
print("\nKey findings:")
print(f"- G82R has {(len(g82r_data) - len(g82_data)) / len(g82_data) * 100:.1f}% more peaks than G82")
print(f"- Mean enrichment increased by {(g82r_data['pValue'].mean() - g82_data['pValue'].mean()) / g82_data['pValue'].mean() * 100:.1f}% in G82R")
print(f"- Peak overlap: ~{overlap} peaks (estimated)")
print(f"- Mean peak width: G82={g82_data['peak_width'].mean():.1f}bp, G82R={g82r_data['peak_width'].mean():.1f}bp")