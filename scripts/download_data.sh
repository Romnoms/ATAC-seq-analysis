#!/bin/bash
# Download ATAC-seq data from GEO

echo "Downloading ATAC-seq peak data from GEO..."
echo "Dataset: GSE130638"
echo "=================================="

# Create data directory
mkdir -p ../data

# Download G82 (parental) peaks
echo "Downloading G82 peaks..."
wget -O ../data/GSM3746025_G82.bed.gz \
  "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3746nnn/GSM3746025/suppl/GSM3746025_G82.bed.gz"

# Download G82R (resistant) peaks  
echo "Downloading G82R peaks..."
wget -O ../data/GSM3746026_G82R.bed.gz \
  "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3746nnn/GSM3746026/suppl/GSM3746026_G82R.bed.gz"

# Download metadata
echo "Downloading metadata..."
wget -O ../data/GSE130638_series_matrix.txt.gz \
  "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE130nnn/GSE130638/matrix/GSE130638_series_matrix.txt.gz"

# Verify downloads
echo ""
echo "Verifying downloads..."
if [ -f "../data/GSM3746025_G82.bed.gz" ] && [ -f "../data/GSM3746026_G82R.bed.gz" ]; then
    echo "✓ Download successful!"
    echo ""
    echo "Files downloaded:"
    ls -lh ../data/*.gz
else
    echo "✗ Download failed. Please check your internet connection and try again."
    exit 1
fi

echo ""
echo "Data download complete!"