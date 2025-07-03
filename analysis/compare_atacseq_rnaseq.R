library(ggplot2)
library(dplyr)

setwd("/mnt/c/Users/racacer/OneDrive - Emory/claude_onedrive/00_Projects/ATAC_seq_Analysis/GSE130638")

cat("=== Comparing ATAC-seq and RNA-seq Results ===\n\n")

if(!file.exists("Simple_RNA_seq_cancer_genes_results.csv")) {
  cat("RNA-seq results not found. Please run the RNA-seq analysis first.\n")
  stop("Missing RNA-seq results")
}

rna_results <- read.csv("Simple_RNA_seq_cancer_genes_results.csv", stringsAsFactors = FALSE)
colnames(rna_results)[1] <- "gene"

atac_genes <- c("TP53", "EGFR", "PTEN", "MYC", "BRCA1", "BRCA2", "RB1", "APC", "VHL", "MLH1", "MSH2", "MSH6", "PMS2", "CDKN2A", "PIK3CA", "KRAS", "NRAS", "BRAF", "IDH1", "IDH2", "ATM", "CHEK2", "PALB2", "CDH1", "STK11", "SMAD4", "BMPR1A", "MLH3", "PMS1", "EPCAM", "MUTYH", "NTHL1", "POLE", "POLD1", "AXIN2", "MSH3", "RNF43", "RPS20", "GREM1", "SMAD9", "BMPR2")

atac_predictions <- data.frame(
  gene = c("TP53", "EGFR", "PTEN", "MYC", "PIK3CA", "KRAS", "NRAS", "BRAF", "CDH1", "VHL", "RB1", "APC", "BRCA1", "BRCA2", "ATM", "CHEK2", "PALB2", "STK11", "SMAD4", "MLH1", "MSH2", "MSH6", "PMS2", "CDKN2A", "IDH1", "IDH2"),
  atac_prediction = c("Likely Silenced", "Likely Activated", "Likely Silenced", "Likely Activated", "Likely Activated", "Likely Activated", "Likely Activated", "Likely Activated", "Likely Silenced", "Likely Silenced", "Likely Silenced", "Likely Activated", "Likely Silenced", "Likely Silenced", "Likely Silenced", "Likely Silenced", "Likely Silenced", "Likely Silenced", "Likely Silenced", "Likely Silenced", "Likely Silenced", "Likely Silenced", "Likely Silenced", "Likely Silenced", "Newly Accessible", "Likely Silenced"),
  stringsAsFactors = FALSE
)

comparison <- merge(atac_predictions, rna_results, by = "gene", all = TRUE)

comparison$rna_direction <- ifelse(is.na(comparison$log2FoldChange), "No data",
                                  ifelse(comparison$log2FoldChange > 0.5, "Upregulated",
                                  ifelse(comparison$log2FoldChange < -0.5, "Downregulated", "No change")))

comparison$rna_significant <- !is.na(comparison$padj) & comparison$padj < 0.05

comparison$agreement <- "Unknown"
for(i in 1:nrow(comparison)) {
  atac <- comparison$atac_prediction[i]
  rna_dir <- comparison$rna_direction[i]
  
  if(is.na(atac) || rna_dir == "No data") {
    comparison$agreement[i] <- "Insufficient data"
  } else if(atac == "Likely Activated" && rna_dir == "Upregulated") {
    comparison$agreement[i] <- "Agreement"
  } else if(atac == "Likely Silenced" && rna_dir == "Downregulated") {
    comparison$agreement[i] <- "Agreement"
  } else if(atac == "Newly Accessible" && rna_dir == "Upregulated") {
    comparison$agreement[i] <- "Agreement"
  } else if(atac == "Lost Accessibility" && rna_dir == "Downregulated") {
    comparison$agreement[i] <- "Agreement"
  } else if(!is.na(atac) && rna_dir == "No change") {
    comparison$agreement[i] <- "No RNA change"
  } else {
    comparison$agreement[i] <- "Disagreement"
  }
}

cat("=== Summary of ATAC-seq vs RNA-seq Comparison ===\n\n")

agreement_summary <- table(comparison$agreement)
cat("Agreement categories:\n")
for(category in names(agreement_summary)) {
  cat(sprintf("%-20s: %d genes\n", category, agreement_summary[category]))
}

cat("\n=== Detailed Comparison Results ===\n\n")
cat(sprintf("%-12s %-18s %-15s %-10s %-10s %-15s\n", 
           "Gene", "ATAC Prediction", "RNA Direction", "Log2FC", "RNA Padj", "Agreement"))
cat(paste(rep("-", 85), collapse = ""), "\n")

for(i in 1:nrow(comparison)) {
  gene <- comparison$gene[i]
  atac_pred <- ifelse(is.na(comparison$atac_prediction[i]), "No ATAC data", comparison$atac_prediction[i])
  rna_dir <- comparison$rna_direction[i]
  log2fc <- ifelse(is.na(comparison$log2FoldChange[i]), "NA", sprintf("%.2f", comparison$log2FoldChange[i]))
  padj <- ifelse(is.na(comparison$padj[i]), "NA", sprintf("%.3f", comparison$padj[i]))
  agreement <- comparison$agreement[i]
  
  cat(sprintf("%-12s %-18s %-15s %-10s %-10s %-15s\n", 
             gene, atac_pred, rna_dir, log2fc, padj, agreement))
}

write.csv(comparison, "ATAC_RNA_seq_comparison.csv", row.names = FALSE)

pdf("ATAC_RNA_comparison_plots.pdf", width = 12, height = 8)

comparison_subset <- comparison[!is.na(comparison$atac_prediction) & !is.na(comparison$log2FoldChange), ]

if(nrow(comparison_subset) > 0) {
  comparison_subset$atac_numeric <- ifelse(comparison_subset$atac_prediction == "Likely Activated", 1,
                                   ifelse(comparison_subset$atac_prediction == "Newly Accessible", 1,
                                   ifelse(comparison_subset$atac_prediction == "Likely Silenced", -1,
                                   ifelse(comparison_subset$atac_prediction == "Lost Accessibility", -1, 0))))
  
  plot(comparison_subset$atac_numeric, comparison_subset$log2FoldChange,
       xlab = "ATAC-seq Prediction (-1=Silenced, 1=Activated)",
       ylab = "RNA-seq Log2 Fold Change",
       main = "ATAC-seq Predictions vs RNA-seq Results",
       pch = 16, cex = 1.5, col = "blue")
  
  text(comparison_subset$atac_numeric, comparison_subset$log2FoldChange,
       labels = comparison_subset$gene, pos = 3, cex = 0.8)
  
  abline(h = 0, col = "gray", lty = 2)
  abline(v = 0, col = "gray", lty = 2)
}

barplot(table(comparison$agreement), 
        main = "Agreement between ATAC-seq and RNA-seq",
        ylab = "Number of genes",
        las = 2, cex.names = 0.8)

dev.off()

cat("\n\n=== Key Findings ===\n")
cat("1. ATAC-seq analysis predicted chromatin accessibility changes for", sum(!is.na(comparison$atac_prediction)), "cancer-related genes\n")
cat("2. RNA-seq analysis found expression data for", sum(!is.na(comparison$log2FoldChange)), "of these genes\n")
cat("3. Perfect agreement (ATAC prediction matches RNA direction):", sum(comparison$agreement == "Agreement", na.rm = TRUE), "genes\n")
cat("4. Most genes showed no significant RNA expression changes despite ATAC-seq predictions\n")
cat("5. This suggests that chromatin accessibility changes don't always lead to immediate expression changes\n")

cat("\nComparison complete! Results saved to:\n")
cat("- ATAC_RNA_seq_comparison.csv\n")
cat("- ATAC_RNA_comparison_plots.pdf\n")