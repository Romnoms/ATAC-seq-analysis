#!/bin/bash
# Run complete ATAC-seq analysis pipeline

echo "ATAC-seq Analysis Pipeline"
echo "=========================="
echo ""

# Set up environment
echo "Setting up environment..."
if command -v conda &> /dev/null; then
    echo "Activating conda environment..."
    conda activate atac-seq 2>/dev/null || echo "Please create conda env: conda env create -f ../environment.yml"
fi

# Check if data exists
if [ ! -f "../data/GSM3746025_G82.bed.gz" ]; then
    echo "Data not found. Running download script..."
    bash download_data.sh
fi

# Run Python analysis
echo ""
echo "Running Python peak analysis..."
cd ../analysis
python plot_atac_peaks.py \
    --g82 ../data/GSM3746025_G82.bed.gz \
    --g82r ../data/GSM3746026_G82R.bed.gz \
    --output ../results/figures

# Run peak comparison
echo ""
echo "Running peak comparison analysis..."
python peak_analysis.py \
    --sample1 ../data/GSM3746025_G82.bed.gz \
    --sample2 ../data/GSM3746026_G82R.bed.gz \
    --name1 G82 \
    --name2 G82R \
    --output ../results/tables/peak_analysis

# Run R analysis
echo ""
echo "Running R analysis notebook..."
Rscript -e "rmarkdown::render('atac_seq_analysis.Rmd', output_dir='../results')"

echo ""
echo "Analysis complete!"
echo "Results saved in:"
echo "  - Figures: ../results/figures/"
echo "  - Tables: ../results/tables/"
echo "  - Report: ../results/atac_seq_analysis.html"