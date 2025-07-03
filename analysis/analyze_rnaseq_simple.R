library(DESeq2)
library(ggplot2)

setwd("/mnt/c/Users/racacer/OneDrive - Emory/claude_onedrive/00_Projects/ATAC_seq_Analysis/GSE130638")

cat("Checking for RNA-seq CSV file...\n")
if (!file.exists("GSE131005_RNASeq.csv")) {
  cat("CSV file not found. Please convert the Excel file to CSV format first.\n")
  cat("You can do this by opening the Excel file and saving it as CSV.\n")
  stop("Missing CSV file")
}

cat("Loading RNA-seq data...\n")
rna_data <- read.csv("GSE131005_RNASeq.csv", stringsAsFactors = FALSE)

cat("Data dimensions:", dim(rna_data), "\n")
cat("Column names:", colnames(rna_data), "\n")

gene_info <- rna_data[, 1:2]
colnames(gene_info) <- c("gene_id", "gene_name")

count_cols <- 3:ncol(rna_data)
count_data <- as.matrix(rna_data[, count_cols])
rownames(count_data) <- gene_info$gene_name

g82_cols <- grep("G82", colnames(count_data), value = TRUE)
count_data <- count_data[, g82_cols]

cat("\nG82 sample columns:", colnames(count_data), "\n")

col_data <- data.frame(
  sample = colnames(count_data),
  condition = ifelse(grepl("D1", colnames(count_data)), "parental", "resistant"),
  row.names = colnames(count_data)
)

cat("\nSample information:\n")
print(col_data)

count_data <- round(count_data)
keep <- rowSums(count_data >= 10) >= 2
count_data_filtered <- count_data[keep, ]
cat("\nGenes remaining after filtering:", nrow(count_data_filtered), "\n")

dds <- DESeqDataSetFromMatrix(
  countData = count_data_filtered,
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
cat("(Comparing G82R resistant vs G82 parental)\n\n")

for(i in 1:nrow(cancer_res)) {
  gene <- rownames(cancer_res)[i]
  log2fc <- cancer_res$log2FoldChange[i]
  padj <- cancer_res$padj[i]
  
  if(!is.na(padj) && padj < 0.05) {
    if(log2fc > 0) {
      cat(sprintf("%-10s: UPREGULATED   (log2FC = %6.2f, padj = %.2e) ***\n", gene, log2fc, padj))
    } else {
      cat(sprintf("%-10s: DOWNREGULATED (log2FC = %6.2f, padj = %.2e) ***\n", gene, log2fc, padj))
    }
  } else {
    cat(sprintf("%-10s: Not significant (log2FC = %6.2f, padj = %.2f)\n", gene, log2fc, ifelse(is.na(padj), 1, padj)))
  }
}

write.csv(res_df, "DESeq2_results_all_genes.csv", row.names = FALSE)
write.csv(sig_genes, "DESeq2_significant_genes.csv", row.names = FALSE)
write.csv(cancer_res, "DESeq2_cancer_genes_results.csv", row.names = TRUE)

pdf("RNA_seq_analysis_plots.pdf", width = 12, height = 10)

plotMA(res, ylim = c(-8, 8), main = "MA Plot: G82R vs G82")

volcano_data <- res_df
volcano_data$significant <- !is.na(volcano_data$padj) & volcano_data$padj < 0.05 & abs(volcano_data$log2FoldChange) > 1
volcano_data$cancer_gene <- volcano_data$gene %in% cancer_genes

p <- ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = significant), alpha = 0.6, size = 1) +
  geom_point(data = volcano_data[volcano_data$cancer_gene & volcano_data$significant, ],
             aes(x = log2FoldChange, y = -log10(pvalue)),
             color = "red", size = 3) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "blue")) +
  theme_minimal() +
  labs(title = "Volcano Plot: G82R vs G82",
       x = "Log2 Fold Change", 
       y = "-Log10 P-value") +
  theme(legend.position = "none")

cancer_sig <- volcano_data[volcano_data$cancer_gene & volcano_data$significant, ]
if(nrow(cancer_sig) > 0) {
  p <- p + geom_text(data = cancer_sig,
                     aes(label = gene), vjust = -0.5, size = 3)
}

print(p)

dev.off()

cat("\nAnalysis complete! Results saved.\n")