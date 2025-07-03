#!/usr/bin/env Rscript

# ATAC-seq Data Visualization for GSE130638
# Comparing G82 (parental) vs G82R (resistant) cells

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(VennDiagram)
library(pheatmap)
library(RColorBrewer)
library(scales)

# Set working directory
setwd("/mnt/c/Users/roman/OneDrive - Emory/claude_onedrive/GSE130638")

# Create output directory for plots
dir.create("plots", showWarnings = FALSE)

# Read the BED files
g82_data <- read.table(gzfile("GSM3746025_G82.bed.gz"), 
                       col.names = c("chr", "start", "end", "name", "score", 
                                    "strand", "signalValue", "pValue", "qValue", 
                                    "peak", "summit", "width", "fold_enrichment"))

g82r_data <- read.table(gzfile("GSM3746026_G82R.bed.gz"), 
                        col.names = c("chr", "start", "end", "name", "score", 
                                     "strand", "signalValue", "pValue", "qValue", 
                                     "peak", "summit", "width", "fold_enrichment"))

# Calculate peak widths
g82_data$peak_width <- g82_data$end - g82_data$start
g82r_data$peak_width <- g82r_data$end - g82r_data$start

# 1. Peak Count Comparison Bar Plot
peak_counts <- data.frame(
  Sample = c("G82\n(Parental)", "G82R\n(Resistant)"),
  Count = c(nrow(g82_data), nrow(g82r_data))
)

p1 <- ggplot(peak_counts, aes(x = Sample, y = Count, fill = Sample)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = format(Count, big.mark = ",")), 
            vjust = -0.5, size = 5, fontface = "bold") +
  scale_fill_manual(values = c("#3498db", "#e74c3c")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(peak_counts$Count) * 1.1),
                     labels = comma) +
  labs(title = "Total ATAC-seq Peaks",
       subtitle = "Comparison between Parental and Resistant Cells",
       y = "Number of Peaks",
       x = "") +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))

ggsave("plots/01_peak_count_comparison.png", p1, width = 8, height = 6, dpi = 300)

# 2. Peak Width Distribution
width_data <- rbind(
  data.frame(Sample = "G82", Width = g82_data$peak_width),
  data.frame(Sample = "G82R", Width = g82r_data$peak_width)
)

p2 <- ggplot(width_data, aes(x = Width, fill = Sample)) +
  geom_histogram(alpha = 0.7, position = "identity", bins = 50) +
  scale_fill_manual(values = c("#3498db", "#e74c3c"),
                    labels = c("G82 (Parental)", "G82R (Resistant)")) +
  scale_x_continuous(limits = c(0, 1500), breaks = seq(0, 1500, 250)) +
  labs(title = "Peak Width Distribution",
       subtitle = "Comparison of ATAC-seq Peak Widths",
       x = "Peak Width (bp)",
       y = "Frequency") +
  theme_minimal() +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        legend.position = "top")

ggsave("plots/02_peak_width_distribution.png", p2, width = 10, height = 6, dpi = 300)

# 3. Chromosome Distribution Heatmap
chr_counts <- rbind(
  g82_data %>% count(chr) %>% mutate(Sample = "G82"),
  g82r_data %>% count(chr) %>% mutate(Sample = "G82R")
) %>%
  filter(chr %in% c(paste0("chr", 1:22), "chrX", "chrY")) %>%
  pivot_wider(names_from = Sample, values_from = n, values_fill = 0)

# Calculate difference and percentage change
chr_counts$Difference <- chr_counts$G82R - chr_counts$G82
chr_counts$Pct_Change <- (chr_counts$Difference / chr_counts$G82) * 100

# Create bar plot for chromosome differences
chr_counts$chr <- factor(chr_counts$chr, 
                        levels = c(paste0("chr", 1:22), "chrX", "chrY"))

p3 <- ggplot(chr_counts, aes(x = chr, y = Pct_Change, fill = Pct_Change > 0)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#e74c3c", "#27ae60"), 
                    labels = c("Decreased", "Increased")) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  labs(title = "Chromosome-specific Peak Changes",
       subtitle = "Percentage change in peak counts (G82R vs G82)",
       x = "Chromosome",
       y = "Percentage Change (%)",
       fill = "Direction") +
  theme_minimal() +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5))

ggsave("plots/03_chromosome_peak_changes.png", p3, width = 10, height = 8, dpi = 300)

# 4. Enrichment Score Distribution
enrich_data <- rbind(
  data.frame(Sample = "G82", Enrichment = g82_data$pValue),
  data.frame(Sample = "G82R", Enrichment = g82r_data$pValue)
)

p4 <- ggplot(enrich_data, aes(x = Sample, y = Enrichment, fill = Sample)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_violin(alpha = 0.5) +
  scale_fill_manual(values = c("#3498db", "#e74c3c")) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(title = "Enrichment Score Distribution",
       subtitle = "MACS2 Peak Enrichment Scores",
       x = "",
       y = "Enrichment Score (-log10 p-value)") +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5))

ggsave("plots/04_enrichment_distribution.png", p4, width = 8, height = 6, dpi = 300)

# 5. Signal Intensity Comparison
signal_data <- rbind(
  data.frame(Sample = "G82", Signal = g82_data$score),
  data.frame(Sample = "G82R", Signal = g82r_data$score)
)

# Create density plot
p5 <- ggplot(signal_data, aes(x = Signal, fill = Sample)) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c("#3498db", "#e74c3c"),
                    labels = c("G82 (Parental)", "G82R (Resistant)")) +
  scale_x_continuous(limits = c(0, 1000), breaks = seq(0, 1000, 200)) +
  labs(title = "Signal Intensity Distribution",
       subtitle = "ATAC-seq Peak Signal Scores",
       x = "Signal Intensity",
       y = "Density") +
  theme_minimal() +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        legend.position = "top")

ggsave("plots/05_signal_intensity_distribution.png", p5, width = 10, height = 6, dpi = 300)

# 6. Create peak overlap analysis
# Convert to GRanges-like format for overlap
g82_regions <- paste0(g82_data$chr, ":", g82_data$start, "-", g82_data$end)
g82r_regions <- paste0(g82r_data$chr, ":", g82r_data$start, "-", g82r_data$end)

# For simplicity, we'll use exact matching (in real analysis, use GenomicRanges)
overlap_count <- length(intersect(g82_regions, g82r_regions))
g82_unique <- length(g82_regions) - overlap_count
g82r_unique <- length(g82r_regions) - overlap_count

# Create Venn diagram
png("plots/06_peak_overlap_venn.png", width = 800, height = 600, res = 150)
venn.plot <- draw.pairwise.venn(
  area1 = length(g82_regions),
  area2 = length(g82r_regions),
  cross.area = overlap_count,
  category = c("G82 (Parental)", "G82R (Resistant)"),
  fill = c("#3498db", "#e74c3c"),
  alpha = 0.5,
  cex = 2,
  cat.cex = 1.5,
  cat.pos = c(-20, 20),
  cat.dist = 0.05
)
dev.off()

# 7. Summary statistics plot
stats_summary <- data.frame(
  Metric = rep(c("Total Peaks", "Mean Width", "Mean Enrichment", "Mean Signal"), each = 2),
  Sample = rep(c("G82", "G82R"), 4),
  Value = c(
    nrow(g82_data), nrow(g82r_data),
    mean(g82_data$peak_width), mean(g82r_data$peak_width),
    mean(g82_data$pValue), mean(g82r_data$pValue),
    mean(g82_data$score), mean(g82r_data$score)
  )
)

# Normalize values for comparison
stats_summary <- stats_summary %>%
  group_by(Metric) %>%
  mutate(Normalized = (Value - min(Value)) / (max(Value) - min(Value))) %>%
  ungroup()

p7 <- ggplot(stats_summary, aes(x = Metric, y = Sample, fill = Normalized)) +
  geom_tile() +
  geom_text(aes(label = round(Value, 1)), size = 5, color = "white", fontface = "bold") +
  scale_fill_gradient2(low = "#3498db", mid = "white", high = "#e74c3c", 
                       midpoint = 0.5, limits = c(0, 1)) +
  labs(title = "Summary Statistics Comparison",
       subtitle = "Key metrics between G82 and G82R samples",
       x = "",
       y = "",
       fill = "Relative\nValue") +
  theme_minimal() +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("plots/07_summary_statistics_heatmap.png", p7, width = 10, height = 4, dpi = 300)

# 8. Top peaks comparison
top_g82 <- g82_data %>% 
  arrange(desc(pValue)) %>% 
  head(20) %>%
  mutate(Sample = "G82", Rank = 1:20)

top_g82r <- g82r_data %>% 
  arrange(desc(pValue)) %>% 
  head(20) %>%
  mutate(Sample = "G82R", Rank = 1:20)

top_peaks <- rbind(top_g82, top_g82r)

p8 <- ggplot(top_peaks, aes(x = Rank, y = pValue, color = Sample)) +
  geom_point(size = 3) +
  geom_line(aes(group = Sample), size = 1) +
  scale_color_manual(values = c("#3498db", "#e74c3c"),
                     labels = c("G82 (Parental)", "G82R (Resistant)")) +
  labs(title = "Top 20 Peaks by Enrichment Score",
       subtitle = "Comparison of highest scoring peaks",
       x = "Peak Rank",
       y = "Enrichment Score (-log10 p-value)") +
  theme_minimal() +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        legend.position = "top")

ggsave("plots/08_top_peaks_comparison.png", p8, width = 10, height = 6, dpi = 300)

# Print summary
cat("\nATAC-seq Analysis Plots Generated Successfully!\n")
cat("===============================================\n")
cat("Output directory: ./plots/\n")
cat("\nGenerated plots:\n")
cat("1. Peak count comparison\n")
cat("2. Peak width distribution\n")
cat("3. Chromosome-specific changes\n")
cat("4. Enrichment score distribution\n")
cat("5. Signal intensity distribution\n")
cat("6. Peak overlap Venn diagram\n")
cat("7. Summary statistics heatmap\n")
cat("8. Top peaks comparison\n")
cat("\nKey findings:\n")
cat(sprintf("- G82R has %.1f%% more peaks than G82\n", 
            (nrow(g82r_data) - nrow(g82_data)) / nrow(g82_data) * 100))
cat(sprintf("- Mean enrichment increased by %.1f%% in G82R\n",
            (mean(g82r_data$pValue) - mean(g82_data$pValue)) / mean(g82_data$pValue) * 100))
cat(sprintf("- Peak overlap: ~%d peaks (estimated)\n", overlap_count))