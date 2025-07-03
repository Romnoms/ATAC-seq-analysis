setwd("/mnt/c/Users/racacer/OneDrive - Emory/claude_onedrive/00_Projects/ATAC_seq_Analysis/GSE130638")

cat("Loading RNA-seq data...\n")
rna_data <- read.csv("GSE131005_RNASeq.csv", stringsAsFactors = FALSE)

cat("Data dimensions:", dim(rna_data), "\n")

g82_cols <- grep("G82", colnames(rna_data))
count_data <- rna_data[, g82_cols]

colnames(count_data) <- gsub(".*=", "", colnames(count_data))
colnames(count_data) <- gsub("-human_align", "", colnames(count_data))

cat("G82 sample columns:", colnames(count_data), "\n")

gene_summary <- rna_data[!duplicated(rna_data$Gene.Symbol), c("Gene.Symbol", "Gene.ID")]
gene_summary <- gene_summary[gene_summary$Gene.Symbol != "", ]

count_summary <- aggregate(count_data, by = list(rna_data$Gene.Symbol), FUN = sum)
count_summary <- count_summary[count_summary$Group.1 != "", ]
rownames(count_summary) <- count_summary$Group.1
count_summary <- count_summary[, -1]

count_matrix <- as.matrix(count_summary)

sample_info <- data.frame(
  sample = colnames(count_matrix),
  condition = ifelse(grepl("D18", colnames(count_matrix)), "resistant", "parental"),
  stringsAsFactors = FALSE
)

cat("\nSample information:\n")
print(sample_info)

parental_samples <- sample_info$sample[sample_info$condition == "parental"]
resistant_samples <- sample_info$sample[sample_info$condition == "resistant"]

cat("Parental samples:", parental_samples, "\n")
cat("Resistant samples:", resistant_samples, "\n")

if(length(parental_samples) == 0 | length(resistant_samples) == 0) {
  stop("Cannot find both parental and resistant samples")
}

keep <- rowSums(count_matrix >= 5) >= 2
count_filtered <- count_matrix[keep, ]
cat("\nGenes remaining after filtering:", nrow(count_filtered), "\n")

parental_counts <- count_filtered[, parental_samples, drop = FALSE]
resistant_counts <- count_filtered[, resistant_samples, drop = FALSE]

parental_mean <- rowMeans(parental_counts + 1)
resistant_mean <- rowMeans(resistant_counts + 1)

log2fc <- log2(resistant_mean / parental_mean)

p_values <- apply(count_filtered, 1, function(x) {
  parental_vals <- x[parental_samples]
  resistant_vals <- x[resistant_samples]
  
  if(all(parental_vals == 0) && all(resistant_vals == 0)) {
    return(1)
  }
  
  tryCatch({
    t.test(resistant_vals, parental_vals)$p.value
  }, error = function(e) {
    return(1)
  })
})

padj <- p.adjust(p_values, method = "BH")

results <- data.frame(
  gene = rownames(count_filtered),
  baseMean = (parental_mean + resistant_mean) / 2,
  log2FoldChange = log2fc,
  pvalue = p_values,
  padj = padj,
  stringsAsFactors = FALSE
)

results <- results[order(results$padj), ]

sig_genes <- results[!is.na(results$padj) & results$padj < 0.05 & abs(results$log2FoldChange) > 1, ]
cat("\nSignificant genes (padj < 0.05, |log2FC| > 1):", nrow(sig_genes), "\n")

if(file.exists("cancer_genes.txt")) {
  cancer_genes <- read.table("cancer_genes.txt", header = TRUE, stringsAsFactors = FALSE)$gene
} else {
  cancer_genes <- c("TP53", "EGFR", "PTEN", "MYC", "BRCA1", "BRCA2", "RB1", "APC", "VHL", "MLH1", "MSH2", "MSH6", "PMS2", "CDKN2A", "PIK3CA", "KRAS", "NRAS", "BRAF", "IDH1", "IDH2", "ATM", "CHEK2", "PALB2", "CDH1", "STK11", "SMAD4", "BMPR1A", "MLH3", "PMS1", "EPCAM", "MUTYH", "NTHL1", "POLE", "POLD1", "AXIN2", "MSH3", "RNF43", "RPS20", "GREM1", "SMAD9", "BMPR2")
}
cancer_genes_found <- intersect(cancer_genes, results$gene)
cancer_res <- results[results$gene %in% cancer_genes_found, ]
cancer_res <- cancer_res[order(cancer_res$log2FoldChange, decreasing = TRUE), ]

cat("\n=== RNA-seq Results for Cancer-Related Genes ===\n")
cat("(Comparing G82 resistant vs parental using simple t-test)\n\n")

for(i in 1:nrow(cancer_res)) {
  gene <- cancer_res$gene[i]
  log2fc <- cancer_res$log2FoldChange[i]
  padj <- cancer_res$padj[i]
  
  if(!is.na(padj) && padj < 0.05 && abs(log2fc) > 1) {
    if(log2fc > 0) {
      cat(sprintf("%-12s: UPREGULATED   (log2FC = %6.2f, padj = %.3f) *** SIGNIFICANT\n", gene, log2fc, padj))
    } else {
      cat(sprintf("%-12s: DOWNREGULATED (log2FC = %6.2f, padj = %.3f) *** SIGNIFICANT\n", gene, log2fc, padj))
    }
  } else {
    status <- "Not significant"
    if(!is.na(padj) && padj < 0.05 && abs(log2fc) <= 1) {
      status <- "Sig. but small FC"
    }
    cat(sprintf("%-12s: %-18s (log2FC = %6.2f, padj = %.3f)\n", gene, status, log2fc, ifelse(is.na(padj), 1, padj)))
  }
}

write.csv(results, "Simple_RNA_seq_results_all_genes.csv", row.names = FALSE)
write.csv(sig_genes, "Simple_RNA_seq_significant_genes.csv", row.names = FALSE)
write.csv(cancer_res, "Simple_RNA_seq_cancer_genes_results.csv", row.names = FALSE)

pdf("Simple_RNA_seq_plots.pdf", width = 12, height = 8)

plot(results$log2FoldChange, -log10(results$pvalue),
     xlab = "Log2 Fold Change", ylab = "-Log10 P-value",
     main = "Volcano Plot: G82 Resistant vs Parental",
     pch = 16, col = "gray", cex = 0.5)

sig_idx <- !is.na(results$padj) & results$padj < 0.05 & abs(results$log2FoldChange) > 1
points(results$log2FoldChange[sig_idx], -log10(results$pvalue[sig_idx]),
       col = "red", pch = 16, cex = 0.8)

cancer_idx <- results$gene %in% cancer_genes & sig_idx
if(sum(cancer_idx) > 0) {
  points(results$log2FoldChange[cancer_idx], -log10(results$pvalue[cancer_idx]),
         col = "blue", pch = 16, cex = 1.2)
  text(results$log2FoldChange[cancer_idx], -log10(results$pvalue[cancer_idx]),
       labels = results$gene[cancer_idx], pos = 3, cex = 0.7)
}

abline(v = c(-1, 1), col = "black", lty = 2)
abline(h = -log10(0.05), col = "black", lty = 2)

dev.off()

cat("\nSimple differential analysis complete!\n")
cat("Results saved to:\n")
cat("- Simple_RNA_seq_results_all_genes.csv\n")
cat("- Simple_RNA_seq_significant_genes.csv\n")
cat("- Simple_RNA_seq_cancer_genes_results.csv\n")
cat("- Simple_RNA_seq_plots.pdf\n")