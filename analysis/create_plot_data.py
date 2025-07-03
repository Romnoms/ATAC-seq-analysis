#!/usr/bin/env python3

import json
import csv

# Read the analysis data
with open('peak_analysis_data.json', 'r') as f:
    data = json.load(f)

# Create CSV files for easy plotting in Excel/R/Python

# 1. Peak counts comparison
with open('plot_data_peak_counts.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Sample', 'Peak_Count'])
    writer.writerow(['G82_Parental', data['peak_counts']['G82']])
    writer.writerow(['G82R_Resistant', data['peak_counts']['G82R']])

# 2. Peak width statistics
with open('plot_data_width_stats.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Statistic', 'G82_Parental', 'G82R_Resistant', 'Percent_Change'])
    for stat in ['mean', 'median', 'q1', 'q3']:
        g82_val = data['width_stats']['G82'][stat]
        g82r_val = data['width_stats']['G82R'][stat]
        pct_change = ((g82r_val - g82_val) / g82_val) * 100
        writer.writerow([stat, g82_val, g82r_val, pct_change])

# 3. Chromosome distribution
with open('plot_data_chromosome_dist.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Chromosome', 'G82_Count', 'G82R_Count', 'Difference', 'Percent_Change'])
    for chr_info in data['chr_distribution']:
        writer.writerow([
            chr_info['chr'], 
            chr_info['g82'], 
            chr_info['g82r'], 
            chr_info['diff'], 
            chr_info['pct_change']
        ])

# 4. Enrichment score statistics
with open('plot_data_enrichment_stats.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Statistic', 'G82_Parental', 'G82R_Resistant', 'Percent_Change'])
    for stat in ['mean', 'median', 'q1', 'q3']:
        g82_val = data['enrichment_stats']['G82'][stat]
        g82r_val = data['enrichment_stats']['G82R'][stat]
        pct_change = ((g82r_val - g82_val) / g82_val) * 100
        writer.writerow([stat, g82_val, g82r_val, pct_change])

# 5. Signal intensity statistics
with open('plot_data_signal_stats.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Statistic', 'G82_Parental', 'G82R_Resistant', 'Percent_Change'])
    for stat in ['mean', 'median', 'q1', 'q3']:
        g82_val = data['signal_stats']['G82'][stat]
        g82r_val = data['signal_stats']['G82R'][stat]
        pct_change = ((g82r_val - g82_val) / g82_val) * 100
        writer.writerow([stat, g82_val, g82r_val, pct_change])

# 6. Overlap analysis
with open('plot_data_overlap.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Category', 'Count', 'Percentage'])
    total_peaks = data['overlap']['total'] + data['overlap']['g82_unique'] + data['overlap']['g82r_unique']
    writer.writerow(['G82_Specific', data['overlap']['g82_unique'], 
                    data['overlap']['g82_unique'] / total_peaks * 100])
    writer.writerow(['G82R_Specific', data['overlap']['g82r_unique'], 
                    data['overlap']['g82r_unique'] / total_peaks * 100])
    writer.writerow(['Overlapping', data['overlap']['total'], 
                    data['overlap']['total'] / total_peaks * 100])

# Create ASCII bar charts for quick visualization
def create_ascii_bar_chart(values, labels, title, width=50):
    print(f"\n{title}")
    print("=" * len(title))
    
    max_val = max(values)
    for i, (label, value) in enumerate(zip(labels, values)):
        bar_length = int((value / max_val) * width)
        bar = "█" * bar_length
        print(f"{label:20} {bar} {value:,.0f}")

# Generate ASCII visualizations
print("ATAC-seq Analysis - Text-based Visualizations")
print("=" * 50)

# Peak counts
create_ascii_bar_chart(
    [data['peak_counts']['G82'], data['peak_counts']['G82R']],
    ['G82 (Parental)', 'G82R (Resistant)'],
    "Peak Counts Comparison"
)

# Chromosome changes (top 10 by absolute change)
chr_data = sorted(data['chr_distribution'], key=lambda x: abs(x['pct_change']), reverse=True)[:10]
create_ascii_bar_chart(
    [abs(x['pct_change']) for x in chr_data],
    [f"{x['chr']} ({x['pct_change']:+.1f}%)" for x in chr_data],
    "Top 10 Chromosome Changes (Absolute %)"
)

# Statistics comparison
stats_comparison = []
labels = []

for metric in ['Peak Count', 'Mean Width', 'Mean Enrichment', 'Mean Signal']:
    if metric == 'Peak Count':
        g82_val = data['peak_counts']['G82']
        g82r_val = data['peak_counts']['G82R']
    elif metric == 'Mean Width':
        g82_val = data['width_stats']['G82']['mean']
        g82r_val = data['width_stats']['G82R']['mean']
    elif metric == 'Mean Enrichment':
        g82_val = data['enrichment_stats']['G82']['mean']
        g82r_val = data['enrichment_stats']['G82R']['mean']
    else:  # Mean Signal
        g82_val = data['signal_stats']['G82']['mean']
        g82r_val = data['signal_stats']['G82R']['mean']
    
    pct_change = ((g82r_val - g82_val) / g82_val) * 100
    stats_comparison.append(pct_change)
    labels.append(f"{metric}")

create_ascii_bar_chart(
    [abs(x) for x in stats_comparison],
    [f"{label} ({val:+.1f}%)" for label, val in zip(labels, stats_comparison)],
    "Percentage Changes in Key Metrics (G82R vs G82)"
)

print(f"\nData files created for external plotting:")
print("- plot_data_peak_counts.csv")
print("- plot_data_width_stats.csv") 
print("- plot_data_chromosome_dist.csv")
print("- plot_data_enrichment_stats.csv")
print("- plot_data_signal_stats.csv")
print("- plot_data_overlap.csv")
print("- ATAC_seq_analysis_report.html (interactive)")

print("\nTo create publication-quality plots:")
print("1. Import CSV files into R, Python, or Excel")
print("2. Open the HTML file in a web browser for interactive plots")
print("3. Use the provided R script template (plot_atac_analysis.R) if R packages are available")