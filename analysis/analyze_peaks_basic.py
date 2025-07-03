#!/usr/bin/env python3

import gzip
import json
from collections import defaultdict, Counter

def read_bed_file(filename):
    """Read BED file and return list of peak data"""
    peaks = []
    with gzip.open(filename, 'rt') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 13:
                peak = {
                    'chr': parts[0],
                    'start': int(parts[1]),
                    'end': int(parts[2]),
                    'name': parts[3],
                    'score': int(parts[4]),
                    'pValue': float(parts[7]),
                    'width': int(parts[2]) - int(parts[1]),
                    'fold_enrichment': float(parts[12])
                }
                peaks.append(peak)
    return peaks

print("Reading ATAC-seq peak files...")
g82_peaks = read_bed_file('GSM3746025_G82.bed.gz')
g82r_peaks = read_bed_file('GSM3746026_G82R.bed.gz')

# 1. Basic statistics
print("\n=== BASIC STATISTICS ===")
print(f"G82 (Parental): {len(g82_peaks):,} peaks")
print(f"G82R (Resistant): {len(g82r_peaks):,} peaks")
print(f"Difference: {len(g82r_peaks) - len(g82_peaks):,} peaks ({(len(g82r_peaks) - len(g82_peaks))/len(g82_peaks)*100:.1f}% increase)")

# 2. Peak width analysis
g82_widths = [p['width'] for p in g82_peaks]
g82r_widths = [p['width'] for p in g82r_peaks]

def calculate_stats(values):
    values.sort()
    n = len(values)
    return {
        'mean': sum(values) / n,
        'median': values[n // 2],
        'q1': values[n // 4],
        'q3': values[3 * n // 4],
        'min': values[0],
        'max': values[-1]
    }

print("\n=== PEAK WIDTH STATISTICS ===")
g82_width_stats = calculate_stats(g82_widths)
g82r_width_stats = calculate_stats(g82r_widths)

print("G82 (Parental):")
for key, value in g82_width_stats.items():
    print(f"  {key}: {value:.1f} bp")

print("\nG82R (Resistant):")
for key, value in g82r_width_stats.items():
    print(f"  {key}: {value:.1f} bp")

# 3. Chromosome distribution
g82_chr_counts = Counter(p['chr'] for p in g82_peaks)
g82r_chr_counts = Counter(p['chr'] for p in g82r_peaks)

print("\n=== CHROMOSOME DISTRIBUTION ===")
print("Chr\tG82\tG82R\tDiff\t%Change")
print("-" * 40)

chr_order = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
chr_data = []

for chr in chr_order:
    if chr in g82_chr_counts or chr in g82r_chr_counts:
        g82_count = g82_chr_counts.get(chr, 0)
        g82r_count = g82r_chr_counts.get(chr, 0)
        diff = g82r_count - g82_count
        pct_change = (diff / g82_count * 100) if g82_count > 0 else 0
        print(f"{chr}\t{g82_count}\t{g82r_count}\t{diff:+d}\t{pct_change:+.1f}%")
        chr_data.append({
            'chr': chr,
            'g82': g82_count,
            'g82r': g82r_count,
            'diff': diff,
            'pct_change': pct_change
        })

# 4. Enrichment score analysis
g82_pvalues = [p['pValue'] for p in g82_peaks]
g82r_pvalues = [p['pValue'] for p in g82r_peaks]

print("\n=== ENRICHMENT SCORE STATISTICS ===")
g82_pvalue_stats = calculate_stats(g82_pvalues)
g82r_pvalue_stats = calculate_stats(g82r_pvalues)

print("G82 (Parental):")
for key, value in g82_pvalue_stats.items():
    print(f"  {key}: {value:.2f}")

print("\nG82R (Resistant):")
for key, value in g82r_pvalue_stats.items():
    print(f"  {key}: {value:.2f}")

# 5. Signal intensity analysis
g82_scores = [p['score'] for p in g82_peaks]
g82r_scores = [p['score'] for p in g82r_peaks]

print("\n=== SIGNAL INTENSITY STATISTICS ===")
g82_score_stats = calculate_stats(g82_scores)
g82r_score_stats = calculate_stats(g82r_scores)

print("G82 (Parental):")
for key, value in g82_score_stats.items():
    print(f"  {key}: {value:.1f}")

print("\nG82R (Resistant):")
for key, value in g82r_score_stats.items():
    print(f"  {key}: {value:.1f}")

# 6. Peak overlap analysis (simplified)
g82_regions = set(f"{p['chr']}:{p['start']}-{p['end']}" for p in g82_peaks)
g82r_regions = set(f"{p['chr']}:{p['start']}-{p['end']}" for p in g82r_peaks)

overlap = len(g82_regions.intersection(g82r_regions))
g82_unique = len(g82_regions) - overlap
g82r_unique = len(g82r_regions) - overlap

print("\n=== PEAK OVERLAP ANALYSIS ===")
print(f"Overlapping peaks: {overlap:,}")
print(f"G82-specific peaks: {g82_unique:,}")
print(f"G82R-specific peaks: {g82r_unique:,}")
print(f"Percentage of G82R peaks overlapping with G82: {overlap/len(g82r_regions)*100:.1f}%")

# 7. Top peaks comparison
print("\n=== TOP 10 PEAKS BY ENRICHMENT ===")
g82_top = sorted(g82_peaks, key=lambda x: x['pValue'], reverse=True)[:10]
g82r_top = sorted(g82r_peaks, key=lambda x: x['pValue'], reverse=True)[:10]

print("\nG82 (Parental) Top 10:")
print("Rank\tChr\tStart\tEnd\tEnrichment")
for i, peak in enumerate(g82_top, 1):
    print(f"{i}\t{peak['chr']}\t{peak['start']}\t{peak['end']}\t{peak['pValue']:.2f}")

print("\nG82R (Resistant) Top 10:")
print("Rank\tChr\tStart\tEnd\tEnrichment")
for i, peak in enumerate(g82r_top, 1):
    print(f"{i}\t{peak['chr']}\t{peak['start']}\t{peak['end']}\t{peak['pValue']:.2f}")

# Save data for external plotting
plot_data = {
    'peak_counts': {'G82': len(g82_peaks), 'G82R': len(g82r_peaks)},
    'width_stats': {'G82': g82_width_stats, 'G82R': g82r_width_stats},
    'chr_distribution': chr_data,
    'enrichment_stats': {'G82': g82_pvalue_stats, 'G82R': g82r_pvalue_stats},
    'signal_stats': {'G82': g82_score_stats, 'G82R': g82r_score_stats},
    'overlap': {'total': overlap, 'g82_unique': g82_unique, 'g82r_unique': g82r_unique}
}

with open('peak_analysis_data.json', 'w') as f:
    json.dump(plot_data, f, indent=2)

print("\n=== ANALYSIS COMPLETE ===")
print("Data saved to peak_analysis_data.json for external plotting")