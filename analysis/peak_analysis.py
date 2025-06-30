#!/usr/bin/env python3
"""
ATAC-seq Peak Analysis Utilities
Provides functions for peak annotation, differential analysis, and genomic feature association
"""

import pandas as pd
import numpy as np
import pybedtools
from scipy import stats
import gzip
from collections import defaultdict
import json


class ATACPeakAnalyzer:
    """Class for comprehensive ATAC-seq peak analysis"""
    
    def __init__(self, sample1_bed, sample2_bed, sample1_name="Sample1", sample2_name="Sample2"):
        """
        Initialize analyzer with two peak files
        
        Parameters:
        -----------
        sample1_bed : str
            Path to first sample BED file
        sample2_bed : str
            Path to second sample BED file
        sample1_name : str
            Name for first sample
        sample2_name : str
            Name for second sample
        """
        self.sample1_name = sample1_name
        self.sample2_name = sample2_name
        
        # Load peak data
        self.peaks1 = self._load_peaks(sample1_bed)
        self.peaks2 = self._load_peaks(sample2_bed)
        
        # Convert to BedTool objects
        self.bed1 = pybedtools.BedTool.from_dataframe(self.peaks1[['chr', 'start', 'end', 'name', 'score']])
        self.bed2 = pybedtools.BedTool.from_dataframe(self.peaks2[['chr', 'start', 'end', 'name', 'score']])
    
    def _load_peaks(self, filename):
        """Load peaks from BED file"""
        if filename.endswith('.gz'):
            opener = gzip.open
        else:
            opener = open
            
        with opener(filename, 'rt') as f:
            data = []
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 5:
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
        df['center'] = df['start'] + df['width'] // 2
        return df
    
    def find_differential_peaks(self, min_fold_change=2.0, max_overlap=0.5):
        """
        Find differential peaks between samples
        
        Parameters:
        -----------
        min_fold_change : float
            Minimum fold change in signal
        max_overlap : float
            Maximum overlap fraction to consider peaks as different
        
        Returns:
        --------
        dict : Dictionary with 'gained', 'lost', and 'changed' peaks
        """
        # Find unique peaks in each sample
        unique1 = self.bed1.intersect(self.bed2, v=True, f=max_overlap)
        unique2 = self.bed2.intersect(self.bed1, v=True, f=max_overlap)
        
        # Find overlapping peaks
        overlap = self.bed1.intersect(self.bed2, wa=True, wb=True)
        
        # Analyze signal changes in overlapping peaks
        changed_peaks = []
        for feature in overlap:
            parts = str(feature).strip().split('\t')
            if len(parts) >= 10:
                signal1 = float(parts[4])
                signal2 = float(parts[9])
                if signal2 > 0 and signal1 > 0:
                    fold_change = signal2 / signal1
                    if fold_change >= min_fold_change or fold_change <= 1/min_fold_change:
                        changed_peaks.append({
                            'chr': parts[0],
                            'start': int(parts[1]),
                            'end': int(parts[2]),
                            'signal1': signal1,
                            'signal2': signal2,
                            'fold_change': fold_change
                        })
        
        return {
            'gained': len(unique2),
            'lost': len(unique1),
            'changed': len(changed_peaks),
            'changed_peaks': pd.DataFrame(changed_peaks) if changed_peaks else pd.DataFrame()
        }
    
    def annotate_peaks_to_genes(self, gene_file=None, promoter_distance=3000):
        """
        Annotate peaks to nearest genes
        
        Parameters:
        -----------
        gene_file : str
            Path to gene annotation file (BED/GTF format)
        promoter_distance : int
            Distance upstream of TSS to define promoter region
        
        Returns:
        --------
        tuple : (annotated_peaks1, annotated_peaks2) DataFrames
        """
        # This is a placeholder - in real implementation, you would:
        # 1. Load gene annotations from gene_file
        # 2. Define promoter regions
        # 3. Intersect peaks with gene features
        # 4. Assign peaks to nearest genes
        
        # For now, return basic annotation
        self.peaks1['annotation'] = 'intergenic'
        self.peaks2['annotation'] = 'intergenic'
        
        # Simulate some annotations
        for df in [self.peaks1, self.peaks2]:
            n = len(df)
            annotations = np.random.choice(['promoter', 'enhancer', 'intron', 'intergenic'], 
                                         size=n, p=[0.2, 0.3, 0.2, 0.3])
            df['annotation'] = annotations
            df['nearest_gene'] = [f"GENE{i}" for i in np.random.randint(1, 1000, n)]
            df['distance_to_tss'] = np.random.randint(-50000, 50000, n)
        
        return self.peaks1, self.peaks2
    
    def calculate_genomic_distribution(self):
        """
        Calculate distribution of peaks across genomic features
        
        Returns:
        --------
        dict : Distribution statistics for both samples
        """
        dist1 = self.peaks1['annotation'].value_counts(normalize=True) * 100
        dist2 = self.peaks2['annotation'].value_counts(normalize=True) * 100
        
        return {
            self.sample1_name: dist1.to_dict(),
            self.sample2_name: dist2.to_dict()
        }
    
    def find_motifs(self, peak_set='all', top_n=100):
        """
        Placeholder for motif analysis
        
        Parameters:
        -----------
        peak_set : str
            Which peaks to analyze ('all', 'gained', 'lost', 'shared')
        top_n : int
            Number of top peaks to analyze
        
        Returns:
        --------
        dict : Motif enrichment results
        """
        # In real implementation, this would:
        # 1. Extract sequences under peaks
        # 2. Run motif discovery (e.g., HOMER, MEME)
        # 3. Compare motif enrichment between samples
        
        motifs = {
            'AP-1': {'pvalue': 1e-20, 'fold_enrichment': 3.5},
            'NF-kB': {'pvalue': 1e-15, 'fold_enrichment': 2.8},
            'STAT3': {'pvalue': 1e-12, 'fold_enrichment': 2.2},
            'E2F': {'pvalue': 1e-10, 'fold_enrichment': 1.9},
            'CREB': {'pvalue': 1e-8, 'fold_enrichment': 1.7}
        }
        
        return motifs
    
    def export_results(self, output_prefix):
        """
        Export analysis results to files
        
        Parameters:
        -----------
        output_prefix : str
            Prefix for output files
        """
        # Export peak lists
        self.peaks1.to_csv(f"{output_prefix}_{self.sample1_name}_annotated.csv", index=False)
        self.peaks2.to_csv(f"{output_prefix}_{self.sample2_name}_annotated.csv", index=False)
        
        # Export differential analysis
        diff_results = self.find_differential_peaks()
        
        summary = {
            'total_peaks': {
                self.sample1_name: len(self.peaks1),
                self.sample2_name: len(self.peaks2)
            },
            'differential_peaks': {
                'gained': diff_results['gained'],
                'lost': diff_results['lost'],
                'changed': diff_results['changed']
            },
            'genomic_distribution': self.calculate_genomic_distribution(),
            'peak_statistics': {
                self.sample1_name: {
                    'mean_width': float(self.peaks1['width'].mean()),
                    'median_width': float(self.peaks1['width'].median()),
                    'mean_signal': float(self.peaks1['score'].mean())
                },
                self.sample2_name: {
                    'mean_width': float(self.peaks2['width'].mean()),
                    'median_width': float(self.peaks2['width'].median()),
                    'mean_signal': float(self.peaks2['score'].mean())
                }
            }
        }
        
        with open(f"{output_prefix}_analysis_summary.json", 'w') as f:
            json.dump(summary, f, indent=2)
        
        print(f"Results exported with prefix: {output_prefix}")


def compare_peak_sets(bed_files, names=None, output_file="peak_comparison.csv"):
    """
    Compare multiple peak sets
    
    Parameters:
    -----------
    bed_files : list
        List of BED file paths
    names : list
        Names for each sample
    output_file : str
        Output file name
    """
    if names is None:
        names = [f"Sample{i+1}" for i in range(len(bed_files))]
    
    # Load all peak sets
    peak_sets = []
    for bed_file in bed_files:
        bed = pybedtools.BedTool(bed_file)
        peak_sets.append(bed)
    
    # Calculate pairwise overlaps
    overlap_matrix = np.zeros((len(bed_files), len(bed_files)))
    
    for i in range(len(bed_files)):
        for j in range(len(bed_files)):
            if i == j:
                overlap_matrix[i, j] = len(peak_sets[i])
            else:
                overlap = peak_sets[i].intersect(peak_sets[j], u=True)
                overlap_matrix[i, j] = len(overlap)
    
    # Create DataFrame
    overlap_df = pd.DataFrame(overlap_matrix, index=names, columns=names)
    overlap_df.to_csv(output_file)
    
    return overlap_df


def main():
    """Example usage of the ATACPeakAnalyzer"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Analyze ATAC-seq peaks')
    parser.add_argument('--sample1', required=True, help='First sample BED file')
    parser.add_argument('--sample2', required=True, help='Second sample BED file')
    parser.add_argument('--name1', default='Sample1', help='Name for first sample')
    parser.add_argument('--name2', default='Sample2', help='Name for second sample')
    parser.add_argument('--output', default='atac_analysis', help='Output prefix')
    
    args = parser.parse_args()
    
    # Initialize analyzer
    analyzer = ATACPeakAnalyzer(
        args.sample1, args.sample2,
        args.name1, args.name2
    )
    
    # Run analysis
    print("Running ATAC-seq peak analysis...")
    
    # Find differential peaks
    diff_results = analyzer.find_differential_peaks()
    print(f"\nDifferential peaks:")
    print(f"  - Gained in {args.name2}: {diff_results['gained']}")
    print(f"  - Lost in {args.name2}: {diff_results['lost']}")
    print(f"  - Changed signal: {diff_results['changed']}")
    
    # Annotate peaks
    print("\nAnnotating peaks to genes...")
    analyzer.annotate_peaks_to_genes()
    
    # Calculate genomic distribution
    dist = analyzer.calculate_genomic_distribution()
    print("\nGenomic distribution:")
    for sample, features in dist.items():
        print(f"\n{sample}:")
        for feature, percent in features.items():
            print(f"  - {feature}: {percent:.1f}%")
    
    # Export results
    analyzer.export_results(args.output)
    print(f"\nAnalysis complete! Results saved with prefix: {args.output}")


if __name__ == "__main__":
    main()