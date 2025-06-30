# ATAC-seq Analysis Pipeline

## Overview
This repository contains a comprehensive ATAC-seq (Assay for Transposase-Accessible Chromatin using sequencing) analysis pipeline for analyzing chromatin accessibility changes between parental and drug-resistant cancer cells.

### Dataset: GSE130638
- **Study**: ATAC-seq on parental (G82) and acquired resistant (G82R) brain tumor cells
- **Treatment**: G82R cells were derived by continuous exposure to CAY10566 inhibitor for ~3 weeks
- **Platform**: Illumina HiSeq 2500
- **Reference Genome**: hg19
- **Publication**: PMID 33568479

## Key Findings

### Global Chromatin Accessibility Changes
- **G82 (Parental)**: 54,414 peaks
- **G82R (Resistant)**: 60,338 peaks (+10.9%)
- **Peak Overlap**: 73.06% of G82R peaks overlap with G82
- **Unique to G82R**: 16,256 new accessible regions

### Peak Characteristics
| Metric | G82 | G82R | Change |
|--------|-----|------|--------|
| Mean Width (bp) | 344.58 | 359.39 | +4.3% |
| Median Width (bp) | 268 | 272 | +1.5% |
| Mean Enrichment | 17.09 | 19.64 | +14.9% |
| Mean Signal | 139.3 | 165.6 | +18.9% |

## Repository Structure

```
ATAC_git/
├── README.md                    # This file
├── requirements.txt             # Python dependencies
├── environment.yml              # Conda environment
├── analysis/
│   ├── atac_seq_analysis.Rmd   # Main R analysis notebook
│   ├── plot_atac_peaks.py      # Python visualization script
│   └── peak_analysis.py        # Peak analysis utilities
├── data/
│   └── README.md               # Data download instructions
├── results/
│   ├── figures/                # Generated plots
│   └── tables/                 # Analysis results
└── scripts/
    ├── download_data.sh        # Data download script
    └── run_analysis.sh         # Full pipeline script
```

## Installation

### Prerequisites
- Python 3.8+
- R 4.0+
- Conda (recommended)

### Setup Python Environment
```bash
# Create conda environment
conda env create -f environment.yml
conda activate atac-seq

# Or use pip
pip install -r requirements.txt
```

### Setup R Environment
```R
# Install required R packages
source("scripts/install_packages.R")
```

## Usage

### Quick Start
```bash
# 1. Download data
bash scripts/download_data.sh

# 2. Run complete analysis
bash scripts/run_analysis.sh

# 3. Generate plots
python analysis/plot_atac_peaks.py
```

### Step-by-Step Analysis

1. **Data Preparation**
   ```bash
   # Download from GEO
   wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3746nnn/GSM3746025/suppl/GSM3746025_G82.bed.gz
   wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3746nnn/GSM3746026/suppl/GSM3746026_G82R.bed.gz
   ```

2. **Run R Analysis**
   ```R
   # Open RStudio and run
   rmarkdown::render("analysis/atac_seq_analysis.Rmd")
   ```

3. **Generate Visualizations**
   ```bash
   python analysis/plot_atac_peaks.py --input data/ --output results/figures/
   ```

## Analysis Pipeline

### 1. Quality Control
- FastQC for raw read quality assessment
- Adapter trimming with TrimGalore
- Alignment statistics evaluation

### 2. Peak Calling
- Alignment with Bowtie2 to hg19
- Duplicate removal with Picard
- Peak calling with MACS2

### 3. Differential Analysis
- Peak overlap analysis
- Signal intensity comparison
- Chromosome-specific changes
- Gene accessibility analysis

### 4. Visualization
- Peak count comparisons
- Peak width distributions
- Chromosome-specific heatmaps
- Enrichment score profiles
- Venn diagrams for peak overlap

## Results

The analysis reveals significant chromatin remodeling in drug-resistant cells:
- Increased overall chromatin accessibility
- Broader peak widths indicating extended accessible regions
- Higher signal intensities suggesting stronger ATAC-seq signal
- Gain of ~16,000 new regulatory regions
- Chromosome-specific responses (Chr3, Chr6, Chr11 most affected)

## Citation

If you use this pipeline, please cite:
- Original study: [Insert publication details]
- This pipeline: [Your repository URL]

## License

MIT License - see LICENSE file for details

## Contact

For questions or issues, please open an issue on GitHub or contact [your email].