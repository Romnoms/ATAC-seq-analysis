#!/usr/bin/env Rscript

# Analyze gene accessibility changes based on ATAC-seq peaks
# This script identifies genes that may be silenced or activated based on chromatin accessibility

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

# Define key genes with their TSS coordinates (hg38)
# These are important cancer-related genes
key_genes <- data.frame(
  gene = c("TP53", "EGFR", "PTEN", "MYC", "CDKN2A", "RB1", "PIK3CA", "BRAF", 
           "KRAS", "ATM", "BRCA1", "BRCA2", "VHL", "APC", "SMAD4", "NF1", 
           "MLH1", "MSH2", "ERBB2", "ALK", "RET", "MET", "FGFR1", "FGFR2",
           "IDH1", "IDH2", "TERT", "ATRX", "H3F3A", "HIST1H3B", "PDGFRA", "PDGFRB",
           "NOTCH1", "PTCH1", "SMO", "CTNNB1", "AKT1", "MTOR", "PBRM1", "BAP1"),
  chr = c("chr17", "chr7", "chr10", "chr8", "chr9", "chr13", "chr3", "chr7",
          "chr12", "chr11", "chr17", "chr13", "chr3", "chr5", "chr18", "chr17",
          "chr3", "chr2", "chr17", "chr2", "chr10", "chr7", "chr8", "chr10",
          "chr2", "chr15", "chr5", "chrX", "chr1", "chr1", "chr4", "chr5",
          "chr9", "chr9", "chr7", "chr3", "chr14", "chr1", "chr3", "chr3"),
  tss = c(7661779, 55019017, 87863438, 127735434, 21968199, 48303747, 179148114, 140719327,
          25227733, 108236969, 43044295, 32314861, 10141776, 112707544, 50190785, 29421945,
          36993455, 47403046, 39688344, 29220830, 43114478, 116966907, 38140658, 123237848,
          208248388, 90017361, 1295162, 77135065, 228041214, 149390804, 54859001, 149351319,
          136546836, 98271958, 129147922, 41001585, 104765224, 11743719, 52302103, 45985687),
  stringsAsFactors = FALSE
)

# Function to find peaks near TSS (within 2kb upstream and 500bp downstream)
find_peaks_near_tss <- function(peaks_df, gene_info, window_upstream=2000, window_downstream=500) {
  nearby_peaks <- data.frame()
  
  for(i in 1:nrow(gene_info)) {
    gene_chr <- gene_info$chr[i]
    gene_tss <- gene_info$tss[i]
    gene_name <- gene_info$gene[i]
    
    # Define promoter region
    promoter_start <- gene_tss - window_upstream
    promoter_end <- gene_tss + window_downstream
    
    # Find peaks in this region
    chr_peaks <- peaks_df[peaks_df$chr == gene_chr, ]
    overlapping <- chr_peaks[chr_peaks$start <= promoter_end & chr_peaks$end >= promoter_start, ]
    
    if(nrow(overlapping) > 0) {
      overlapping$gene <- gene_name
      overlapping$distance_to_tss <- pmin(abs(overlapping$start - gene_tss), 
                                         abs(overlapping$end - gene_tss))
      nearby_peaks <- rbind(nearby_peaks, overlapping)
    }
  }
  
  return(nearby_peaks)
}

# Find peaks near TSS for both samples
cat("Finding peaks near gene promoters...\n")
g82_promoter_peaks <- find_peaks_near_tss(g82_df, key_genes)
g82r_promoter_peaks <- find_peaks_near_tss(g82r_df, key_genes)

# Summarize accessibility by gene
g82_gene_summary <- aggregate(cbind(score, pValue) ~ gene, data=g82_promoter_peaks, FUN=max)
colnames(g82_gene_summary)[2:3] <- c("g82_max_score", "g82_max_enrichment")

g82r_gene_summary <- aggregate(cbind(score, pValue) ~ gene, data=g82r_promoter_peaks, FUN=max)
colnames(g82r_gene_summary)[2:3] <- c("g82r_max_score", "g82r_max_enrichment")

# Merge summaries
gene_comparison <- merge(g82_gene_summary, g82r_gene_summary, by="gene", all=TRUE)
gene_comparison[is.na(gene_comparison)] <- 0

# Calculate changes
gene_comparison$score_change <- gene_comparison$g82r_max_score - gene_comparison$g82_max_score
gene_comparison$enrichment_change <- gene_comparison$g82r_max_enrichment - gene_comparison$g82_max_enrichment
gene_comparison$score_fold_change <- (gene_comparison$g82r_max_score + 1) / (gene_comparison$g82_max_score + 1)

# Categorize genes
gene_comparison$category <- "No Change"
gene_comparison$category[gene_comparison$score_change > 50 | gene_comparison$enrichment_change > 5] <- "Likely Activated"
gene_comparison$category[gene_comparison$score_change < -50 | gene_comparison$enrichment_change < -5] <- "Likely Silenced"
gene_comparison$category[gene_comparison$g82_max_score == 0 & gene_comparison$g82r_max_score > 0] <- "Newly Accessible"
gene_comparison$category[gene_comparison$g82_max_score > 0 & gene_comparison$g82r_max_score == 0] <- "Lost Accessibility"

# Sort by absolute change
gene_comparison <- gene_comparison[order(abs(gene_comparison$score_change), decreasing=TRUE), ]

# Create separate visualizations for each chart
par(family="sans")

# 1. Bar chart - Top Genes by Accessibility Change
png("atac_plots/gene_barplot.png", width=1920, height=1200, res=150)
par(mar=c(8,5,3,1.6))

# Ensure we have at least 20 genes to plot
genes_to_plot <- min(20, nrow(gene_comparison))
genes_subset <- gene_comparison[1:genes_to_plot, ]

bp <- barplot(genes_subset$score_change, 
              names.arg=genes_subset$gene, 
              las=2, 
              col=ifelse(genes_subset$score_change > 0, "#e74c3c", "#3498db"),
              main="Top Genes by Accessibility Change",
              ylab="",
              cex.names=1.12,
              cex.main=1.6,
              cex.lab=0.9,
              cex.axis=1.28)

# Add y-axis title manually with reduced size and left positioning
mtext("Signal Score Change (G82R - G82)", side=2, line=5, cex=0.9)

abline(h=0, lty=2, lwd=2.4)

# Add value labels on bars with proportionally smaller font
text(bp, genes_subset$score_change + sign(genes_subset$score_change) * max(abs(genes_subset$score_change)) * 0.05,
     round(genes_subset$score_change, 0), cex=0.96, adj=0.5, font=2)

dev.off()

# 2. Pie chart - Gene Category Distribution
png("atac_plots/gene_piechart.png", width=1200, height=1200, res=150)
par(mar=c(8,8,8,8))

category_counts <- table(gene_comparison$category)
# Create better colors and labels
colors <- c("Likely Activated"="#e74c3c", "Likely Silenced"="#3498db", 
            "No Change"="lightgray", "Newly Accessible"="#f39c12", 
            "Lost Accessibility"="#9b59b6")
pie_colors <- colors[names(category_counts)]

pie(category_counts, 
    col=pie_colors,
    main="Gene Category Distribution",
    cex.main=2.4,
    labels=paste(names(category_counts), "\n(", category_counts, ")", sep=""),
    cex=2.0)

dev.off()

# 3. Scatter plot - Promoter Accessibility Comparison
png("atac_plots/gene_scatterplot.png", width=1600, height=2000, res=150)
par(mar=c(10,10,8,4))

# Increase x and y axis by 25%
max_x <- max(gene_comparison$g82_max_score, na.rm=TRUE)
max_y <- max(gene_comparison$g82r_max_score, na.rm=TRUE)
xlim_scatter <- c(0, max_x * 1.25)
ylim_scatter <- c(0, max_y * 1.25)

plot(gene_comparison$g82_max_score, gene_comparison$g82r_max_score,
     pch=19, col=rgb(0,0,0,0.7), cex=3.0,
     xlab="G82 Max Signal Score", ylab="G82R Max Signal Score",
     main="Promoter Accessibility Comparison",
     cex.main=3.0, cex.lab=2.5, cex.axis=2.2,
     xlim=xlim_scatter, ylim=ylim_scatter)
abline(a=0, b=1, lty=2, col="red", lwd=6)

# Add labels for top changed genes only, with much larger fonts
top_changed <- head(gene_comparison[order(abs(gene_comparison$score_change), decreasing=TRUE), ], 8)
for(i in 1:nrow(top_changed)) {
  x <- top_changed$g82_max_score[i]
  y <- top_changed$g82r_max_score[i]
  label <- top_changed$gene[i]
  
  # Special positioning for CTNNB1 - place to the left
  if(label == "CTNNB1") {
    text(x, y, label, cex=2.0, pos=2, offset=1.0, col="black", font=2)
  } else {
    # Position other labels to avoid overlap
    text(x, y, label, cex=2.0, pos=4, offset=1.0, col="black", font=2)
  }
}

dev.off()

# 4. Heatmap - Signal Intensity
png("atac_plots/gene_heatmap.png", width=1000, height=1200, res=150)
par(mar=c(8,12,6,12))

top_genes <- head(gene_comparison[gene_comparison$category != "No Change", ], 10)  # Reduced further for readability
if(nrow(top_genes) > 0) {
  gene_matrix <- as.matrix(top_genes[, c("g82_max_score", "g82r_max_score")])
  rownames(gene_matrix) <- top_genes$gene
  colnames(gene_matrix) <- c("G82", "G82R")
  
  # Use actual values instead of normalized for clearer interpretation
  max_val <- max(gene_matrix)
  
  image(1:2, 1:nrow(gene_matrix), t(gene_matrix), 
        col=colorRampPalette(c("white", "yellow", "orange", "red"))(100),
        xlab="", ylab="", axes=FALSE,
        main="Signal Intensity Heatmap",
        cex.main=2.0)
  
  # Add axes with appropriate fonts
  axis(1, at=1:2, labels=c("G82", "G82R"), cex.axis=1.8)
  axis(2, at=1:nrow(gene_matrix), labels=rownames(gene_matrix), las=2, cex.axis=1.4)
  
  # Add value labels in cells with appropriate font
  for(i in 1:nrow(gene_matrix)) {
    for(j in 1:ncol(gene_matrix)) {
      text(j, i, round(gene_matrix[i,j], 0), cex=1.2, col="black", font=2)
    }
  }
  
  # Add color scale legend positioned below the heatmap
  legend("bottom", legend=c("High", "Medium", "Low"), 
         fill=c("red", "orange", "white"), cex=1.2, horiz=TRUE, 
         inset=c(0, -0.15), xpd=TRUE)
}

dev.off()

# Create a separate detailed bar chart showing all genes with changes
png("atac_plots/detailed_gene_changes.png", width=2000, height=1000, res=150)
par(mar=c(14,8,6,3), family="sans")

# Filter genes with actual changes
changed_genes <- gene_comparison[gene_comparison$category != "No Change", ]
changed_genes <- changed_genes[order(changed_genes$score_change, decreasing=TRUE), ]

if(nrow(changed_genes) > 0) {
  bp <- barplot(changed_genes$score_change,
                names.arg=changed_genes$gene,
                las=2,
                col=ifelse(changed_genes$score_change > 0, "#e74c3c", "#3498db"),
                main="All Genes with Significant Accessibility Changes",
                ylab="",
                cex.names=1.8,
                cex.main=2.2,
                cex.lab=1.0,
                cex.axis=1.8)
  
  # Add y-axis title manually with reduced size and left positioning
  mtext("Signal Score Change (G82R - G82)", side=2, line=6, cex=1.0)
  
  abline(h=0, lty=2, lwd=4)
  
  # Add category labels with larger font
  text(bp, changed_genes$score_change + sign(changed_genes$score_change) * max(abs(changed_genes$score_change)) * 0.07,
       substr(changed_genes$category, 1, 1), cex=1.4, font=2)
  
  # Add legend with larger font
  legend("topright", 
         legend=c("Likely Activated", "Newly Accessible", "Likely Silenced", "Lost Accessibility"),
         fill=c("#e74c3c", "#f39c12", "#3498db", "#9b59b6"),
         cex=1.4)
}

dev.off()

# Save detailed results
write.csv(gene_comparison, "atac_plots/gene_accessibility_summary.csv", row.names=FALSE)

# Create summary for HTML
activated_genes <- gene_comparison[gene_comparison$category == "Likely Activated", "gene"]
silenced_genes <- gene_comparison[gene_comparison$category == "Likely Silenced", "gene"]
newly_accessible <- gene_comparison[gene_comparison$category == "Newly Accessible", "gene"]
lost_accessibility <- gene_comparison[gene_comparison$category == "Lost Accessibility", "gene"]

cat("\n=== Gene Accessibility Analysis Summary ===\n")
cat("Likely Activated Genes:", length(activated_genes), "\n")
if(length(activated_genes) > 0) cat("  ", paste(head(activated_genes, 10), collapse=", "), "\n")
cat("\nLikely Silenced Genes:", length(silenced_genes), "\n")
if(length(silenced_genes) > 0) cat("  ", paste(head(silenced_genes, 10), collapse=", "), "\n")
cat("\nNewly Accessible Genes:", length(newly_accessible), "\n")
if(length(newly_accessible) > 0) cat("  ", paste(head(newly_accessible, 10), collapse=", "), "\n")
cat("\nLost Accessibility:", length(lost_accessibility), "\n")
if(length(lost_accessibility) > 0) cat("  ", paste(head(lost_accessibility, 10), collapse=", "), "\n")

# Save summary for HTML integration
summary_list <- list(
  activated = activated_genes,
  silenced = silenced_genes,
  newly_accessible = newly_accessible,
  lost_accessibility = lost_accessibility,
  total_analyzed = nrow(key_genes),
  genes_with_peaks = nrow(gene_comparison[gene_comparison$g82_max_score > 0 | gene_comparison$g82r_max_score > 0, ])
)

if("jsonlite" %in% installed.packages()[,"Package"]) {
  jsonlite::write_json(summary_list, "atac_plots/gene_accessibility_summary.json", pretty=TRUE)
}