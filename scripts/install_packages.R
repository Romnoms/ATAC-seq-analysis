#!/usr/bin/env Rscript
# Install required R packages for ATAC-seq analysis

cat("Installing R packages for ATAC-seq analysis...\n")
cat("============================================\n\n")

# Function to install packages
install_if_missing <- function(pkg, bioc = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing:", pkg, "\n")
    if (bioc) {
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    } else {
      install.packages(pkg, repos = "https://cran.rstudio.com/")
    }
  } else {
    cat("Already installed:", pkg, "\n")
  }
}

# Install BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cran.rstudio.com/")
}

# CRAN packages
cat("\nInstalling CRAN packages...\n")
cran_packages <- c(
  "tidyverse",
  "data.table",
  "ggplot2",
  "pheatmap",
  "RColorBrewer",
  "knitr",
  "rmarkdown",
  "DT",
  "plotly",
  "ggforce"
)

for (pkg in cran_packages) {
  install_if_missing(pkg, bioc = FALSE)
}

# Bioconductor packages
cat("\nInstalling Bioconductor packages...\n")
bioc_packages <- c(
  "GenomicRanges",
  "rtracklayer",
  "ChIPseeker",
  "TxDb.Hsapiens.UCSC.hg19.knownGene",
  "org.Hs.eg.db",
  "clusterProfiler",
  "DESeq2",
  "edgeR",
  "ComplexHeatmap",
  "EnhancedVolcano"
)

for (pkg in bioc_packages) {
  install_if_missing(pkg, bioc = TRUE)
}

cat("\nPackage installation complete!\n")
cat("You can now run the analysis with: Rscript -e \"rmarkdown::render('analysis/atac_seq_analysis.Rmd')\"\n")