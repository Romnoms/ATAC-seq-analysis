#!/usr/bin/env Rscript

# Create ATAC-seq average plots and heatmaps from BED peak files
# This script creates pseudo-signal visualizations using peak data

library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(scales)
library(gridExtra)

# Function to read BED files
read_bed_file <- function(filename) {
  # Read the gzipped BED file
  df <- read.table(gzfile(filename), sep="\t", stringsAsFactors=FALSE)
  
  # Assign column names based on MACS2 BED format
  colnames(df) <- c('chr', 'start', 'end', 'name', 'score', 'strand', 
                    'signalValue', 'pValue', 'qValue', 'peak', 'summit', 
                    'width', 'fold_enrichment')[1:ncol(df)]
  
  return(df)
}

# Create chromosome heatmap
create_chromosome_heatmap <- function(g82_df, g82r_df) {
  # Calculate peak counts and mean enrichment by chromosome
  g82_chr <- g82_df %>%
    group_by(chr) %>%
    summarise(
      peak_count = n(),
      mean_enrichment = mean(pValue),
      mean_signal = mean(score)
    ) %>%
    mutate(sample = "G82")
  
  g82r_chr <- g82r_df %>%
    group_by(chr) %>%
    summarise(
      peak_count = n(),
      mean_enrichment = mean(pValue),
      mean_signal = mean(score)
    ) %>%
    mutate(sample = "G82R")
  
  # Combine data
  chr_data <- rbind(g82_chr, g82r_chr)
  
  # Create matrix for heatmap
  peak_matrix <- chr_data %>%
    select(chr, sample, peak_count) %>%
    pivot_wider(names_from = sample, values_from = peak_count) %>%
    column_to_rownames("chr")
  
  enrichment_matrix <- chr_data %>%
    select(chr, sample, mean_enrichment) %>%
    pivot_wider(names_from = sample, values_from = mean_enrichment) %>%
    column_to_rownames("chr")
  
  # Sort chromosomes
  chr_order <- c(paste0("chr", 1:22), "chrX", "chrY")
  peak_matrix <- peak_matrix[chr_order[chr_order %in% rownames(peak_matrix)], ]
  enrichment_matrix <- enrichment_matrix[chr_order[chr_order %in% rownames(enrichment_matrix)], ]
  
  # Create heatmaps
  png("atac_chromosome_heatmap.png", width=1400, height=1000, res=100)
  par(mfrow=c(1,2))
  
  # Peak count heatmap
  pheatmap(t(peak_matrix), 
           main="Peak Counts by Chromosome",
           cluster_rows=FALSE, 
           cluster_cols=FALSE,
           display_numbers=TRUE,
           number_format="%.0f",
           color=colorRampPalette(c("white", "yellow", "orange", "red"))(100))
  
  # Calculate percentage changes
  pct_change <- (peak_matrix[,"G82R"] - peak_matrix[,"G82"]) / peak_matrix[,"G82"] * 100
  pct_change_matrix <- as.matrix(pct_change)
  colnames(pct_change_matrix) <- "Peak Count Change (%)"
  
  pheatmap(t(pct_change_matrix),
           main="Percentage Change in Peak Count (G82R vs G82)",
           cluster_rows=FALSE, 
           cluster_cols=FALSE,
           display_numbers=TRUE,
           number_format="%.1f",
           color=colorRampPalette(c("blue", "white", "red"))(100),
           breaks=seq(-30, 30, length.out=101))
  
  dev.off()
}

# Create enrichment profile plot
create_enrichment_profile <- function(g82_df, g82r_df, top_n=10000) {
  # Sort by enrichment score
  g82_sorted <- g82_df[order(g82_df$pValue, decreasing=TRUE), ]
  g82r_sorted <- g82r_df[order(g82r_df$pValue, decreasing=TRUE), ]
  
  # Take top peaks
  n_peaks <- min(top_n, nrow(g82_sorted), nrow(g82r_sorted))
  
  # Create rank data
  rank_data <- data.frame(
    rank = 1:n_peaks,
    G82 = g82_sorted$pValue[1:n_peaks],
    G82R = g82r_sorted$pValue[1:n_peaks]
  )
  
  # Calculate moving average
  window <- 50
  rank_data$G82_ma <- zoo::rollmean(rank_data$G82, k=window, fill=NA, align="center")
  rank_data$G82R_ma <- zoo::rollmean(rank_data$G82R, k=window, fill=NA, align="center")
  
  # Create plot
  p1 <- ggplot(rank_data, aes(x=rank)) +
    geom_line(aes(y=G82_ma, color="G82 (Parental)"), size=1.5, na.rm=TRUE) +
    geom_line(aes(y=G82R_ma, color="G82R (Resistant)"), size=1.5, na.rm=TRUE) +
    scale_color_manual(values=c("G82 (Parental)"="#3498db", "G82R (Resistant)"="#e74c3c")) +
    labs(x="Peak Rank", 
         y="Enrichment Score (-log10 p-value)",
         title=paste0("Average Enrichment Profile (Top ", n_peaks, " peaks)"),
         color="Sample") +
    theme_minimal() +
    theme(legend.position="top")
  
  # Peak width vs enrichment scatter plot
  p2 <- ggplot() +
    geom_point(data=g82_df[g82_df$width < 2000 & g82_df$pValue < 50, ], 
               aes(x=width, y=pValue, color="G82"), alpha=0.3, size=0.5) +
    geom_point(data=g82r_df[g82r_df$width < 2000 & g82r_df$pValue < 50, ], 
               aes(x=width, y=pValue, color="G82R"), alpha=0.3, size=0.5) +
    scale_color_manual(values=c("G82"="#3498db", "G82R"="#e74c3c")) +
    labs(x="Peak Width (bp)", 
         y="Enrichment Score (-log10 p-value)",
         title="Peak Width vs Enrichment Score",
         color="Sample") +
    theme_minimal() +
    theme(legend.position="top")
  
  # Save plot
  png("atac_enrichment_profile.png", width=1200, height=1000, res=100)
  grid.arrange(p1, p2, nrow=2)
  dev.off()
}

# Create signal intensity heatmap
create_signal_heatmap <- function(g82_df, g82r_df, n_regions=100) {
  # Get top peaks by signal from each sample
  g82_top <- g82_df[order(g82_df$score, decreasing=TRUE)[1:n_regions], ]
  g82r_top <- g82r_df[order(g82r_df$score, decreasing=TRUE)[1:n_regions], ]
  
  # Create matrix for heatmap
  # Combine unique peaks
  all_peaks <- rbind(
    data.frame(id=paste(g82_top$chr, g82_top$start, g82_top$end, sep="_"),
               G82=g82_top$score,
               G82R=0),
    data.frame(id=paste(g82r_top$chr, g82r_top$start, g82r_top$end, sep="_"),
               G82=0,
               G82R=g82r_top$score)
  )
  
  # Aggregate by peak ID (in case of duplicates)
  peak_matrix <- all_peaks %>%
    group_by(id) %>%
    summarise(G82=max(G82), G82R=max(G82R)) %>%
    column_to_rownames("id") %>%
    as.matrix()
  
  # Transpose for samples as rows
  peak_matrix <- t(peak_matrix)
  
  # Create heatmap
  png("atac_signal_heatmap.png", width=1000, height=600, res=100)
  
  # Scale the data for better visualization
  pheatmap(peak_matrix,
           main=paste0("Signal Intensity Heatmap (Top ", n_regions, " peaks from each sample)"),
           cluster_rows=FALSE,
           cluster_cols=TRUE,
           show_colnames=FALSE,
           color=colorRampPalette(c("white", "yellow", "orange", "red", "darkred"))(100),
           scale="column",
           fontsize=10)
  
  dev.off()
}

# Create genomic distribution plot
create_genomic_distribution <- function(g82_df, g82r_df) {
  # Chromosome lengths (approximate, hg38)
  chr_lengths <- c(
    chr1=248956422, chr2=242193529, chr3=198295559, chr4=190214555,
    chr5=181538259, chr6=170805979, chr7=159345973, chr8=145138636,
    chr9=138394717, chr10=133797422, chr11=135086622, chr12=133275309,
    chr13=114364328, chr14=107043718, chr15=101991189, chr16=90338345,
    chr17=83257441, chr18=80373285, chr19=58617616, chr20=64444167,
    chr21=46709983, chr22=50818468, chrX=156040895, chrY=57227415
  )
  
  # Calculate peak density
  g82_density <- g82_df %>%
    group_by(chr) %>%
    summarise(count = n()) %>%
    mutate(
      length = chr_lengths[chr],
      density = count / length * 1e6,
      sample = "G82"
    )
  
  g82r_density <- g82r_df %>%
    group_by(chr) %>%
    summarise(count = n()) %>%
    mutate(
      length = chr_lengths[chr],
      density = count / length * 1e6,
      sample = "G82R"
    )
  
  density_data <- rbind(g82_density, g82r_density)
  
  # Order chromosomes
  chr_order <- c(paste0("chr", 1:22), "chrX", "chrY")
  density_data$chr <- factor(density_data$chr, levels=chr_order)
  
  # Create density plot
  p1 <- ggplot(density_data, aes(x=chr, y=density, fill=sample)) +
    geom_bar(stat="identity", position="dodge") +
    scale_fill_manual(values=c("G82"="#3498db", "G82R"="#e74c3c")) +
    labs(x="Chromosome", y="Peak Density (peaks per Mb)", 
         title="Peak Density by Chromosome", fill="Sample") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle=45, hjust=1))
  
  # Peak width distribution
  width_data <- rbind(
    data.frame(width=g82_df$width[g82_df$width < 2000], sample="G82"),
    data.frame(width=g82r_df$width[g82r_df$width < 2000], sample="G82R")
  )
  
  p2 <- ggplot(width_data, aes(x=width, fill=sample)) +
    geom_histogram(alpha=0.7, position="identity", bins=50) +
    scale_fill_manual(values=c("G82"="#3498db", "G82R"="#e74c3c")) +
    labs(x="Peak Width (bp)", y="Count", 
         title="Peak Width Distribution", fill="Sample") +
    theme_minimal()
  
  # Save plots
  png("atac_genomic_distribution.png", width=1400, height=600, res=100)
  grid.arrange(p1, p2, ncol=2)
  dev.off()
}

# Main function
main <- function() {
  # Read BED files
  cat("Reading BED files...\n")
  g82_df <- read_bed_file("GSM3746025_G82.bed.gz")
  g82r_df <- read_bed_file("GSM3746026_G82R.bed.gz")
  
  cat(paste0("G82 peaks: ", nrow(g82_df), "\n"))
  cat(paste0("G82R peaks: ", nrow(g82r_df), "\n"))
  
  # Create visualizations
  cat("Creating chromosome heatmap...\n")
  create_chromosome_heatmap(g82_df, g82r_df)
  
  cat("Creating enrichment profile...\n")
  create_enrichment_profile(g82_df, g82r_df)
  
  cat("Creating signal heatmap...\n")
  create_signal_heatmap(g82_df, g82r_df)
  
  cat("Creating genomic distribution plot...\n")
  create_genomic_distribution(g82_df, g82r_df)
  
  cat("All visualizations created successfully!\n")
  
  # Save plot paths for HTML integration
  plot_paths <- list(
    chromosome_heatmap = "atac_chromosome_heatmap.png",
    enrichment_profile = "atac_enrichment_profile.png",
    signal_heatmap = "atac_signal_heatmap.png",
    genomic_distribution = "atac_genomic_distribution.png"
  )
  
  jsonlite::write_json(plot_paths, "atac_plot_paths.json", pretty=TRUE)
}

# Check if required libraries are installed
required_packages <- c("ggplot2", "dplyr", "tidyr", "pheatmap", "RColorBrewer", 
                       "scales", "gridExtra", "zoo", "jsonlite")

missing_packages <- required_packages[!required_packages %in% installed.packages()[,"Package"]]

if(length(missing_packages) > 0) {
  cat("Installing missing packages:", paste(missing_packages, collapse=", "), "\n")
  install.packages(missing_packages, repos="https://cloud.r-project.org/")
}

# Load libraries
suppressMessages({
  library(zoo)
  library(jsonlite)
})

# Run main function
main()