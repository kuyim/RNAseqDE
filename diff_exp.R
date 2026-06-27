# Differential expression: dexamethasone-treated vs untreated airway smooth
# muscle cells (GSE52778), blocking on cell line.
# Input: gene-level count matrix from nf-core/rnaseq (Salmon pseudo-alignment)

library(DESeq2)
library(tidyverse)
library(pheatmap)

# Loads length scaled counts matrix into table
counts_raw <- read.delim(
  "results/airway_rnaseq/salmon/salmon.merged.gene_counts_length_scaled.tsv",
  row.names = 1,
  check.names = FALSE
)

# Allow gene names to be accessed by row name
print(head(counts_raw))
gene_names <- counts_raw[["gene_name"]]
names(gene_names) <- rownames(counts_raw)
counts <- counts_raw %>% dplyr::select(-dplyr::any_of("gene_name"))

# Round counts b/c DESeq2 requires integers
counts <- round(as.matrix(counts))
print(colnames(counts)) # Check colnames before making df

# Create sample metadata table
coldata <- data.frame(
  sample = c("SRX384345", "SRX384346", "SRX384349", "SRX384350",
             "SRX384353", "SRX384354", "SRX384357", "SRX384358"),
  cell_line = c("N61311", "N61311", "N052611", "N052611",
                "N080611", "N080611", "N061011", "N061011"),
  condition = c("untreated", "dexamethasone", "untreated", "dexamethasone",
                "untreated", "dexamethasone", "untreated", "dexamethasone"),
  stringsAsFactors = FALSE
)
rownames(coldata) <- coldata$sample

# Ensure matrix columns are aligned to metadata rows for DESeq2
counts <- counts[, coldata$sample]
stopifnot(all(colnames(counts) == rownames(coldata)))

# Set untreated as reference
coldata$condition <- factor(coldata$condition, levels = c("untreated", "dexamethasone"))
coldata$cell_line <- factor(coldata$cell_line)

# Building DESeqDataSet: needs counts, metadata, and groups
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ cell_line + condition
)

# Filter weakly expressed genes
keep <- rowSums(counts(dds) >= 10) >= 4
dds <- dds[keep, ]
message(sprintf("Num genes retained after filtering: %d", nrow(dds)))


# Run the DESeq
dds <- DESeq(dds)
resultsNames(dds)

dir.create("results/deseq2", showWarnings=FALSE, recursive = TRUE)
vsd <- vst(dds, blind = FALSE) # Variance-stabilizing transform for vis

# Create PCA plots
p <- plotPCA(vsd, intgroup = c("condition", "cell_line")) +
  ggplot2::theme_bw()
ggplot2::ggsave("results/deseq2/pca.png", p, width=7, height=5, dpi=150)


# Sample-distance heatmap: replicates/conditions should cluster
samp_dist <- dist(t(assay(vsd)))
samp_dist_mat <- as.matrix(samp_dist)

# Build "CellLine_condition" labels in the same order as the matrix
srx_ids <- colnames(samp_dist_mat)
cond_short <- ifelse(coldata[srx_ids, "condition"] == "untreated", "unt", "trt")
labels <- paste(coldata[srx_ids, "cell_line"], cond_short, sep = "_")

rownames(samp_dist_mat) <- labels
colnames(samp_dist_mat) <- labels

pheatmap(
  samp_dist_mat,
  clustering_distance_rows = samp_dist,
  clustering_distance_cols = samp_dist,
  filename = "results/deseq2/sample_distances.png"
)


# Extract results
res <- results(dds, contrast = c("condition", "dexamethasone", "untreated"))
summary(res)

# Count and report significantly up/down regulated genes at padj < 0.05
n_tested <- sum(!is.na(res$padj))                       # genes with a valid padj
sig      <- !is.na(res$padj) & res$padj < 0.05
n_sig    <- sum(sig)
n_up     <- sum(sig & res$log2FoldChange > 0)
n_down   <- sum(sig & res$log2FoldChange < 0)
message(sprintf("Significant genes (padj < 0.05): %d (%.1f%% of %d tested)",
                n_sig, 100 * n_sig / n_tested, n_tested))
message(sprintf("  Up   (log2FC > 0): %d (%.1f%%)", n_up,   100 * n_up   / n_tested))
message(sprintf("  Down (log2FC < 0): %d (%.1f%%)", n_down, 100 * n_down / n_tested))

# Shrink log2 fold changes to make ranking/vis more reliable
resLFC <- lfcShrink(dds, coef = "condition_dexamethasone_vs_untreated",
                    type="apeglm")

# MA plot: shrunk log2FC vs mean expression; sig genes (padj < 0.05) highlighted
png("results/deseq2/ma_plot.png", width = 7, height = 5, units = "in", res = 150)
DESeq2::plotMA(resLFC, alpha = 0.05,
               main = "MA plot: dexamethasone vs untreated")
dev.off()

# Assemble, annotate with gene names, order by sig, and save
res_df <- as.data.frame(resLFC) %>%
  tibble::rownames_to_column("gene_id") %>%
  dplyr::mutate(gene_name = gene_names[gene_id]) %>%
  dplyr::relocate(gene_name, .after=gene_id) %>%
  dplyr::arrange(padj)

readr::write_csv(res_df, "results/deseq2/deseq2_dexamethasone_vs_untreated.csv")
message("Done. Results written to csv.")


# Heatmap of the top N Differentially Expressed Genes (DEGs) (by padj), using VST-transformed expression
topN <- 50
top_ids <- res_df %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::slice_head(n = topN) %>%
  dplyr::pull(gene_id)

# VST matrix for the top genes; center each gene so relative changes are visible
top_mat <- assay(vsd)[top_ids, , drop = FALSE]
top_mat <- top_mat - rowMeans(top_mat)


# Label rows with gene names (fall back to gene_id) and columns with cell_line_cond
rownames(top_mat) <- ifelse(is.na(gene_names[top_ids]) | gene_names[top_ids] == "",
                            top_ids, gene_names[top_ids])
colnames(top_mat) <- labels

annotation_col <- data.frame(
  condition = coldata[srx_ids, "condition"],
  cell_line = coldata[srx_ids, "cell_line"],
  row.names = labels
)

pheatmap(
  top_mat,
  annotation_col = annotation_col,
  main = sprintf("Top %d DEGs (VST, centered)", topN),
  filename = "results/deseq2/top_degs_heatmap.png"
)
message(sprintf("Saved MA plot and top-%d DEG heatmap to results/deseq2/.", topN))





