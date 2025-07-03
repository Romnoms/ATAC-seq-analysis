library(DESeq2)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(tidyr)

setwd("/mnt/c/Users/racacer/OneDrive - Emory/claude_onedrive/00_Projects/ATAC_seq_Analysis/GSE130638")

cat("Attempting to load RNA-seq data...\n")

if (!require("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx", repos = "http://cran.r-project.org")
  library(openxlsx)
}

cat("Loading RNA-seq data from Excel file...\n")
rna_data <- read.xlsx("GSE131005_RNASeq.xlsx", sheet = 1)

cat("Data dimensions:", dim(rna_data), "\n")
cat("Column names:", colnames(rna_data), "\n")

gene_info <- rna_data[, 1:2]
colnames(gene_info) <- c("gene_id", "gene_name")

count_data <- rna_data[, 3:ncol(rna_data)]
count_data <- as.matrix(count_data)
rownames(count_data) <- gene_info$gene_name

count_data <- count_data[, grepl("G82", colnames(count_data))]

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

cat("\n=== RNA-seq Results for 40 Cancer-Related Genes ===\n")
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

sig_genes_for_heatmap <- head(sig_genes$gene, 50)
vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)[sig_genes_for_heatmap, ]
mat <- mat - rowMeans(mat)

pheatmap(mat, 
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "euclidean",
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         annotation_col = col_data[, "condition", drop = FALSE],
         main = "Top 50 Differentially Expressed Genes",
         fontsize_row = 8)

cancer_genes_sig <- cancer_res[!is.na(cancer_res$padj) & cancer_res$padj < 0.05, ]
if(nrow(cancer_genes_sig) > 0) {
  cancer_mat <- assay(vsd)[rownames(cancer_genes_sig), ]
  cancer_mat <- cancer_mat - rowMeans(cancer_mat)
  
  pheatmap(cancer_mat,
           clustering_distance_rows = "correlation",
           clustering_distance_cols = "euclidean", 
           cluster_cols = TRUE,
           show_rownames = TRUE,
           show_colnames = TRUE,
           annotation_col = col_data[, "condition", drop = FALSE],
           main = "Significantly Changed Cancer-Related Genes",
           fontsize_row = 10)
}

volcano_data <- res_df %>%
  mutate(
    significant = !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1,
    cancer_gene = gene %in% cancer_genes
  )

ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = significant), alpha = 0.6, size = 1) +
  geom_point(data = volcano_data[volcano_data$cancer_gene & volcano_data$significant, ],
             aes(x = log2FoldChange, y = -log10(pvalue)),
             color = "red", size = 3) +
  geom_text(data = volcano_data[volcano_data$cancer_gene & volcano_data$significant, ],
            aes(label = gene), vjust = -0.5, size = 3) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "blue")) +
  theme_minimal() +
  labs(title = "Volcano Plot: G82R vs G82",
       x = "Log2 Fold Change",
       y = "-Log10 P-value") +
  theme(legend.position = "none")

dev.off()

cat("\nAnalysis complete! Results saved to:\n")
cat("- DESeq2_results_all_genes.csv\n")
cat("- DESeq2_significant_genes.csv\n") 
cat("- DESeq2_cancer_genes_results.csv\n")
cat("- RNA_seq_analysis_plots.pdf\n")