# Data Directory

This directory contains ATAC-seq peak data files. The data files are not included in the repository due to their size.

## Required Data Files

To run the analysis, you need the following files:
- `GSM3746025_G82.bed.gz` - ATAC-seq peaks for G82 (parental) cells
- `GSM3746026_G82R.bed.gz` - ATAC-seq peaks for G82R (resistant) cells

## Download Instructions

### Option 1: Automatic Download
Run the download script from the scripts directory:
```bash
cd scripts
bash download_data.sh
```

### Option 2: Manual Download
Download the files directly from GEO:
1. Visit https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130638
2. Download the supplementary files:
   - GSM3746025_G82.bed.gz
   - GSM3746026_G82R.bed.gz
3. Place them in this directory

## Data Format

The BED files are in MACS2 narrowPeak format with the following columns:
1. **chr** - Chromosome
2. **start** - Peak start position
3. **end** - Peak end position
4. **name** - Peak name
5. **score** - Peak score (0-1000)
6. **strand** - Strand (. for ATAC-seq)
7. **signalValue** - Signal value
8. **pValue** - -log10(p-value)
9. **qValue** - -log10(q-value)
10. **peak** - Peak summit position relative to start

## Data Source

- **GEO Accession**: GSE130638
- **SRA Study**: SRP194611
- **Platform**: Illumina HiSeq 2500
- **Organism**: Homo sapiens
- **Reference**: hg19