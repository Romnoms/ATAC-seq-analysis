#!/usr/bin/env Rscript

# Simple ATAC-seq visualization using base R
# Creates average plots and heatmaps from BED peak files

# Read BED files
g82_df <- read.table(gzfile("GSM3746025_G82.bed.gz"), sep="\t", stringsAsFactors=FALSE)
g82r_df <- read.table(gzfile("GSM3746026_G82R.bed.gz"), sep="\t", stringsAsFactors=FALSE)

# Assign column names
colnames(g82_df) <- c('chr', 'start', 'end', 'name', 'score', 'strand', 
                      'signalValue', 'pValue', 'qValue', 'peak', 'summit', 
                      'width', 'fold_enrichment')[1:ncol(g82_df)]
colnames(g82r_df) <- c('chr', 'start', 'end', 'name', 'score', 'strand', 
                       'signalValue', 'pValue', 'qValue', 'peak', 'summit', 
                       'width', 'fold_enrichment')[1:ncol(g82r_df)]

cat("G82 peaks:", nrow(g82_df), "\n")
cat("G82R peaks:", nrow(g82r_df), "\n")

# Create output directory for plots
dir.create("atac_plots", showWarnings = FALSE)

# 1. Create enrichment profile plot
png("atac_plots/enrichment_profile.png", width=1200, height=900, res=150)
par(mfrow=c(2,1), mar=c(5,5,4,3), family="sans")

# Sort by enrichment and plot
g82_sorted <- sort(g82_df$pValue, decreasing=TRUE)
g82r_sorted <- sort(g82r_df$pValue, decreasing=TRUE)

n_peaks <- min(5000, length(g82_sorted), length(g82r_sorted))
# Increase y-axis by 100% for enrichment score
max_enrichment <- max(c(g82_sorted[1:n_peaks], g82r_sorted[1:n_peaks]))
ylim_enrichment <- c(0, max_enrichment * 2)

plot(1:n_peaks, g82_sorted[1:n_peaks], type='l', col='#3498db', lwd=3,
     xlab="Peak Rank", ylab="Enrichment Score (-log10 p-value)",
     main="Enrichment Score Profile (Top 5000 peaks)",
     cex.main=1.3, cex.lab=1.2, cex.axis=1.1,
     ylim=ylim_enrichment)
lines(1:n_peaks, g82r_sorted[1:n_peaks], col='#e74c3c', lwd=3)
legend("topright", c("G82 (Parental)", "G82R (Resistant)"), 
       col=c('#3498db', '#e74c3c'), lwd=3, cex=1.1)

# Peak width distribution
hist(g82_df$width[g82_df$width < 2000], breaks=50, col=rgb(52/255,152/255,219/255,0.5),
     xlim=c(0,2000), ylim=c(0,8000), main="Peak Width Distribution",
     xlab="Peak Width (bp)", ylab="Frequency",
     cex.main=1.3, cex.lab=1.2, cex.axis=1.1)
hist(g82r_df$width[g82r_df$width < 2000], breaks=50, col=rgb(231/255,76/255,60/255,0.5), add=TRUE)
legend("topright", c("G82", "G82R"), fill=c(rgb(52/255,152/255,219/255,0.5), rgb(231/255,76/255,60/255,0.5)), cex=1.1)
dev.off()

# 2. Create chromosome-based heatmap
png("atac_plots/chromosome_heatmap.png", width=1200, height=800, res=150)
par(family="sans")

# Calculate peaks per chromosome
chr_order <- c(paste0("chr", 1:22), "chrX", "chrY")
g82_chr_counts <- table(factor(g82_df$chr, levels=chr_order))
g82r_chr_counts <- table(factor(g82r_df$chr, levels=chr_order))

# Create matrix for heatmap
chr_matrix <- rbind(g82_chr_counts, g82r_chr_counts)
rownames(chr_matrix) <- c("G82", "G82R")

# Plot heatmap with better margins and sizing
heatmap(as.matrix(chr_matrix), Colv=NA, Rowv=NA, scale="none",
        col=colorRampPalette(c("white", "yellow", "orange", "red"))(100),
        main="Peak Counts by Chromosome", margins=c(6,6),
        cexRow=1.2, cexCol=1.0)
dev.off()

# 3. Create signal intensity heatmap for top peaks
png("atac_plots/signal_heatmap.png", width=1200, height=600, res=150)
par(family="sans")

# Get top 100 peaks by signal
g82_top <- g82_df[order(g82_df$score, decreasing=TRUE)[1:100], ]
g82r_top <- g82r_df[order(g82r_df$score, decreasing=TRUE)[1:100], ]

# Create matrix
signal_matrix <- matrix(0, nrow=2, ncol=200)
signal_matrix[1, 1:100] <- g82_top$score
signal_matrix[2, 101:200] <- g82r_top$score
rownames(signal_matrix) <- c("G82", "G82R")

# Plot heatmap
image(1:200, 1:2, t(signal_matrix), col=colorRampPalette(c("white", "yellow", "red"))(100),
      xlab="Top 100 peaks from each sample", ylab="", yaxt='n',
      main="Signal Intensity Heatmap")
axis(2, at=1:2, labels=c("G82", "G82R"), las=2)
dev.off()

# 4. Create average plot around peak summits
png("atac_plots/average_plot.png", width=1200, height=900, res=150)
par(mfrow=c(2,1), mar=c(5,5,4,3), family="sans")

# Simulate average profile using peak widths and scores
# This is a simplified representation
window_size <- 1000
x <- seq(-window_size, window_size, 10)

# Create gaussian-like profiles based on peak widths
g82_profile <- numeric(length(x))
g82r_profile <- numeric(length(x))

# Use top 1000 peaks
n_top <- min(1000, nrow(g82_df), nrow(g82r_df))
g82_top_peaks <- g82_df[order(g82_df$pValue, decreasing=TRUE)[1:n_top], ]
g82r_top_peaks <- g82r_df[order(g82r_df$pValue, decreasing=TRUE)[1:n_top], ]

# Create average profiles
for(i in 1:length(x)) {
  g82_profile[i] <- mean(exp(-0.5 * (x[i]^2) / (g82_top_peaks$width/4)^2) * g82_top_peaks$score)
  g82r_profile[i] <- mean(exp(-0.5 * (x[i]^2) / (g82r_top_peaks$width/4)^2) * g82r_top_peaks$score)
}

# Increase y-axis by 25% for average signal profile
max_signal <- max(c(g82_profile, g82r_profile))
ylim_signal <- c(0, max_signal * 1.25)

plot(x, g82_profile, type='l', col='#3498db', lwd=3,
     xlab="Distance from Peak Center (bp)", ylab="Average Signal",
     main="Average Signal Profile at Peak Centers",
     cex.main=1.3, cex.lab=1.2, cex.axis=1.1,
     ylim=ylim_signal)
lines(x, g82r_profile, col='#e74c3c', lwd=3)
legend("topright", c("G82", "G82R"), col=c('#3498db', '#e74c3c'), lwd=3, cex=1.1)
abline(v=0, lty=2, col="gray", lwd=2)

# Difference plot with appropriate y-axis scaling
diff_signal <- g82r_profile - g82_profile
max_diff <- max(abs(diff_signal))
ylim_diff <- c(-max_diff * 1.25, max_diff * 1.25)

plot(x, diff_signal, type='l', col='darkgreen', lwd=3,
     xlab="Distance from Peak Center (bp)", ylab="Signal Difference",
     main="Signal Difference (G82R - G82)",
     cex.main=1.3, cex.lab=1.2, cex.axis=1.1,
     ylim=ylim_diff)
abline(h=0, lty=2, col="gray", lwd=2)
abline(v=0, lty=2, col="gray", lwd=2)
dev.off()

cat("All plots created successfully in atac_plots/ directory\n")

# Save summary statistics
summary_stats <- data.frame(
  Metric = c("Total Peaks", "Mean Peak Width", "Mean Enrichment Score", "Mean Signal"),
  G82 = c(nrow(g82_df), mean(g82_df$width), mean(g82_df$pValue), mean(g82_df$score)),
  G82R = c(nrow(g82r_df), mean(g82r_df$width), mean(g82r_df$pValue), mean(g82r_df$score))
)

write.csv(summary_stats, "atac_plots/summary_stats.csv", row.names=FALSE)
cat("Summary statistics saved to atac_plots/summary_stats.csv\n")