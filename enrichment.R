# Pathway / GO enrichment on the dexamethasone-vs-untreated DE results.
# Input : results/deseq2/deseq2_dexamethasone_vs_untreated.csv  (from diff_exp.R)
# Output: results/enrichment/  (tables + plots)

library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)     # human gene annotation database
library(enrichplot)

# Load DE results, use gene_name as the symbol column
res_df <- readr::read_csv("results/deseq2/deseq2_dexamethasone_vs_untreated.csv")
print(head(res_df$gene_name))


# Map gene symbols to ENTREZID, see how many map and how many don't
res_df <- res_df %>% filter(!is.na(gene_name), gene_name != "")

id_map <- bitr(unique(res_df$gene_name),
               fromType = "SYMBOL",
               toType = "ENTREZID",
               OrgDb = org.Hs.eg.db)

res_mapped <- res_df %>%
  inner_join(id_map, by = c("gene_name" = "SYMBOL")) %>%
  filter(!is.na(ENTREZID)) %>%
  distinct(ENTREZID, .keep_all = TRUE)
message(sprintf("Genes mapped to Entrez: %d of %d", nrow(res_mapped), nrow(res_df)))


# Define gene sets: Universe is all tested genes with valid padj, not whole genome
universe <- res_mapped %>% filter(!is.na(padj)) %>% pull(ENTREZID)
sig_genes <- res_mapped %>% filter(!is.na(padj), padj < 0.05) %>% pull(ENTREZID)
message(sprintf("Signficant genes (padj < 0.05) mapped: %d", length(sig_genes)))

dir.create("results/enrichment", showWarnings = FALSE, recursive = TRUE)

save_table <- function(obj, path) {
  df <- as.data.frame(obj)
  if (nrow(df) > 0) readr::write_csv(df, path)
  invisible(df)
}


# Overrepresentation Analysis (ORA): Gene Ontology (GO)
# Checks if DEGs are overrepresented from any particular pathway
ego <- enrichGO(
  gene          = sig_genes,
  universe      = universe,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.10,
  readable      = TRUE
)

ego_df <- save_table(ego, "results/enrichment/go_bp_ora.csv")
if (nrow(ego_df) > 0) {
  p <- dotplot(ego, showCategory = 20) + ggtitle("GO BP over-representation (DEGs)")
  ggsave("results/enrichment/go_bp_ora_dotplot.png", p, width = 9, height = 8, dpi = 150)
}


# ORA: Kyoto Encyclopedia of Genes and Genomes (KEGG) pathways
ekegg <- tryCatch(
  enrichKEGG(gene = sig_genes, universe = universe,
             organism = "hsa", pvalueCutoff = 0.05),
  error = function(e) { message("KEGG ORA skipped: ", conditionMessage(e)); NULL}
)

if (!is.null(ekegg)) {
  ek_df <- save_table(ekegg, "results/enrichment/kegg_ora.csv")
  if (nrow(ek_df) > 0) {
    p <- dotplot(ekegg, showCategory = 20) + ggtitle("KEGG over-representation (DEGs)")
    ggsave("results/enrichment/kegg_ora_dotplot.png", p, width = 9, height = 8, dpi = 150)
  }
}



# Gene Set Enrichment Analysis (GSEA): GO biological process
ranked <- res_mapped %>% filter(!is.na(log2FoldChange)) %>%
          arrange(desc(log2FoldChange))
gene_list <- ranked$log2FoldChange
names(gene_list) <- ranked$ENTREZID
gene_list <- sort(gene_list, decreasing=TRUE)

set.seed(42)   # GSEA uses a stochastic permutation method
gse <- gseGO(
  geneList     = gene_list,
  OrgDb        = org.Hs.eg.db,
  ont          = "BP",
  keyType      = "ENTREZID",
  pvalueCutoff = 0.05,
  seed         = TRUE,
  verbose      = FALSE
)

gse_df <- save_table(gse, "results/enrichment/go_bp_gsea.csv")
if (nrow(gse_df) > 0) {
  # Dotplot split by activated vs suppressed pathways
  p <- dotplot(gse, showCategory = 15, split = ".sign") +
    facet_grid(. ~ .sign) + ggtitle("GO BP GSEA")
  ggsave("results/enrichment/go_bp_gsea_dotplot.png", p, width = 11, height = 9, dpi = 150)
  # Ridgeplot of core-enrichment fold-change distributions
  pr <- ridgeplot(gse, showCategory = 20) + ggtitle("GO BP GSEA (fold-change ridges)")
  ggsave("results/enrichment/go_bp_gsea_ridge.png", pr, width = 10, height = 9, dpi = 150)
}

message("Enrichment complete. Tables and plots written to results/enrichment/.")


# Check where ORA:GO and GSEA:GO agree; KEGG uses a different pathway vocab
# Pull each result into a common schema (NA NES for the ORA methods).
tidy_enrich <- function(obj, analysis, has_nes = FALSE) {
  if (is.null(obj)) return(NULL)
  df <- as.data.frame(obj)
  if (nrow(df) == 0) return(NULL)
  tibble::tibble(
    analysis    = analysis,
    ID          = df$ID,
    Description = df$Description,
    p.adjust    = df$p.adjust,
    NES         = if (has_nes) df$NES else NA_real_,
    direction   = if (has_nes) ifelse(df$NES > 0, "up", "down") else "enriched"
  )
}

summary_all <- dplyr::bind_rows(
  tidy_enrich(ego,   "GO_ORA"),
  tidy_enrich(ekegg, "KEGG_ORA"),
  tidy_enrich(gse,   "GSEA_GO", has_nes = TRUE)
)

if (!is.null(summary_all) && nrow(summary_all) > 0) {
  # Full combined listing of every enriched term from all three analyses
  readr::write_csv(summary_all, "results/enrichment/enrichment_summary_all.csv")
  
  # Headline convergence: GO terms flagged by BOTH GO-ORA and GSEA-GO
  convergence_go <- dplyr::inner_join(
    summary_all %>% dplyr::filter(analysis == "GO_ORA") %>%
      dplyr::select(ID, Description, ORA_padjust = p.adjust),
    summary_all %>% dplyr::filter(analysis == "GSEA_GO") %>%
      dplyr::select(ID, GSEA_padjust = p.adjust, GSEA_NES = NES, direction),
    by = "ID"
  ) %>% dplyr::arrange(ORA_padjust)
  
  if (nrow(convergence_go) > 0) {
    readr::write_csv(convergence_go,
                     "results/enrichment/convergence_go_ora_vs_gsea.csv")
    message(sprintf("GO terms corroborated by BOTH ORA and GSEA: %d",
                    nrow(convergence_go)))
    print(utils::head(convergence_go, 15))
  } else {
    message("No GO terms were shared between ORA and GSEA at current cutoffs.")
  }
} else {
  message("No enriched terms to summarize.")
}

message("Enrichment complete. Tables and plots written to results/enrichment/.")








