if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "http://cran.r-project.org")

if (!require("DESeq2", quietly = TRUE))
  BiocManager::install("DESeq2")

library(DESeq2)

setwd("/mnt/c/Users/racacer/OneDrive - Emory/claude_onedrive/00_Projects/ATAC_seq_Analysis/GSE130638")

cat("Loading RNA-seq data...\n")
rna_data <- read.csv("GSE131005_RNASeq.csv", stringsAsFactors = FALSE)

cat("Data dimensions:", dim(rna_data), "\n")

gene_info <- rna_data[, 1:4]

g82_columns <- c(
  "run1768_lane2_indexD704.D508.G82D18A.human_align",
  "run1768_lane2_indexD705.D507.G82D18B.human_align", 
  "run1768_lane2_indexD708.D503.G82DOB.human_align",
  "run1773_lane92_indexD707.D504.G82DOA.human_align"
)

available_cols <- colnames(rna_data)
cat("Available columns containing G82:\n")
for(col in available_cols) {
  if(grepl("G82", col)) {
    cat(" -", col, "\n")
  }
}

count_cols <- grep("G82", colnames(rna_data))
count_data <- as.matrix(rna_data[, count_cols])

colnames(count_data) <- gsub(".*=", "", colnames(count_data))
colnames(count_data) <- gsub("-human_align", "", colnames(count_data))

cat("\nProcessed column names:", colnames(count_data), "\n")

gene_summary <- rna_data[!duplicated(rna_data$Gene.Symbol), c("Gene.Symbol", "Gene.ID")]
gene_summary <- gene_summary[gene_summary$Gene.Symbol != "", ]

count_summary <- aggregate(count_data, by = list(rna_data$Gene.Symbol), FUN = sum)
count_summary <- count_summary[count_summary$Group.1 != "", ]
rownames(count_summary) <- count_summary$Group.1
count_summary <- count_summary[, -1]

count_matrix <- as.matrix(count_summary)
count_matrix <- round(count_matrix)

col_data <- data.frame(
  sample = colnames(count_matrix),
  condition = ifelse(grepl("D18", colnames(count_matrix)), "resistant", "parental"),
  row.names = colnames(count_matrix)
)

cat("\nSample information:\n")
print(col_data)

keep <- rowSums(count_matrix >= 10) >= 2
count_filtered <- count_matrix[keep, ]
cat("\nGenes remaining after filtering:", nrow(count_filtered), "\n")

dds <- DESeqDataSetFromMatrix(
  countData = count_filtered,
  colData = col_data,
  design = ~ condition
)

dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "resistant", "parental"))

res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df <- res_df[order(res_df$padj), ]

sig_genes <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, ]
cat("\nSignificant genes (padj < 0.05, |log2FC| > 1):", nrow(sig_genes), "\n")

cancer_genes <- read.table("cancer_genes.txt", header = TRUE, stringsAsFactors = FALSE)$gene
cancer_genes_found <- intersect(cancer_genes, rownames(res_df))
cancer_res <- res_df[cancer_genes_found, ]
cancer_res <- cancer_res[order(cancer_res$log2FoldChange, decreasing = TRUE), ]

cat("\n=== RNA-seq Results for Cancer-Related Genes ===\n")
cat("(Comparing G82 resistant vs G82 parental)\n\n")

for(i in 1:nrow(cancer_res)) {
  gene <- rownames(cancer_res)[i]
  log2fc <- cancer_res$log2FoldChange[i]
  padj <- cancer_res$padj[i]
  
  if(!is.na(padj) && padj < 0.05) {
    if(log2fc > 0) {
      cat(sprintf("%-10s: UPREGULATED   (log2FC = %6.2f, padj = %.2e) *** SIGNIFICANT\n", gene, log2fc, padj))
    } else {
      cat(sprintf("%-10s: DOWNREGULATED (log2FC = %6.2f, padj = %.2e) *** SIGNIFICANT\n", gene, log2fc, padj))
    }
  } else {
    cat(sprintf("%-10s: Not significant (log2FC = %6.2f, padj = %.2f)\n", gene, log2fc, ifelse(is.na(padj), 1, padj)))
  }
}

write.csv(res_df, "RNA_seq_DESeq2_results_all_genes.csv", row.names = FALSE)
write.csv(sig_genes, "RNA_seq_DESeq2_significant_genes.csv", row.names = FALSE)
write.csv(cancer_res, "RNA_seq_DESeq2_cancer_genes_results.csv", row.names = TRUE)

cat("\nSaving plots...\n")
pdf("RNA_seq_analysis_plots.pdf", width = 12, height = 8)

plotMA(res, ylim = c(-8, 8), main = "MA Plot: G82 Resistant vs Parental")

dev.off()

cat("\nAnalysis complete!\n")
cat("Results saved to:\n")
cat("- RNA_seq_DESeq2_results_all_genes.csv\n")
cat("- RNA_seq_DESeq2_significant_genes.csv\n")
cat("- RNA_seq_DESeq2_cancer_genes_results.csv\n")
cat("- RNA_seq_analysis_plots.pdf\n")