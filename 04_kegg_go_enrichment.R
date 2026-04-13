# ==============================================================================
# Script: 04_kegg_go_enrichment.R
# Project: Temporal Gene Expression Dynamics of NEC in Human Intestinal
#          Organoid Model (Ellaithi, Gattu Linga et al.)
#
# Purpose:
#   Perform KEGG pathway and Gene Ontology (GO) over-representation analysis
#   on significant DEGs identified by LIMMA (Figures 3-6 in the manuscript).
#
# Method:
#   - Auto-detects gene symbol, logFC, and adjusted P-value columns
#   - Converts HGNC gene symbols to Entrez IDs using org.Hs.eg.db (human)
#   - KEGG pathway enrichment via enrichKEGG (organism: hsa)
#   - GO enrichment via enrichGO across BP, CC, and MF ontologies
#   - BH-adjusted P < 0.05 for enrichment significance
#   - COVID-19-related KEGG terms removed for biological relevance
#   - Generates dot plots and bar plots with wrapped pathway labels
#
# Input:  results/LIMMA_DEG_significant_NEC_vs_Control_{24h,72h}.csv
# Output:
#   - results/KEGG_results_{24h,72h}.csv
#   - results/GO_{BP,CC,MF}_results_{24h,72h}.csv
#   - results/SYMBOL_to_ENTREZID_{24h,72h}.csv
#   - results/KEGG_{dotplot,barplot}_{24h,72h}.{pdf,jpg}
#   - results/GO_{BP,CC,MF}_dotplot_{24h,72h}.{pdf,jpg}
#
# Usage:
#   Set `timepoint` below to "24h" or "72h", then source this script.
#
# References:
#   clusterProfiler: Wu et al., The Innovation (2021)
#   org.Hs.eg.db: Carlson M, Bioconductor annotation package
# ==============================================================================

# ---- Configuration -----------------------------------------------------------
timepoint    <- "24h"  # Change to "72h" for 72-hour analysis
padj_cutoff  <- 0.05
logfc_cutoff <- 1
show_kegg    <- 20
show_go      <- 20
dpi_out      <- 900
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
})

# Define timepoint-specific parameters
if (timepoint == "24h") {
  deg_file   <- "results/LIMMA_DEG_significant_NEC_vs_Control_24h.csv"
  out_suffix <- "24h"
} else if (timepoint == "72h") {
  deg_file   <- "results/LIMMA_DEG_significant_NEC_vs_Control_72h.csv"
  out_suffix <- "72h"
} else {
  stop("timepoint must be '24h' or '72h'")
}

# ---- Helper functions --------------------------------------------------------
pick_col <- function(df, candidates) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

apply_readable_theme <- function(p, wrap_width = 55, y_text_size = 11) {
  p +
    scale_y_discrete(labels = function(x) str_wrap(x, width = wrap_width)) +
    theme_bw() +
    theme(
      axis.text.y  = element_text(size = y_text_size, face = "bold"),
      axis.text.x  = element_text(size = 11),
      axis.title.x = element_text(size = 12, face = "bold"),
      axis.title.y = element_text(size = 12, face = "bold"),
      plot.title   = element_text(face = "bold", size = 14, hjust = 0.5),
      legend.title = element_text(size = 11, face = "bold"),
      legend.text  = element_text(size = 11),
      plot.margin  = margin(t = 10, r = 10, b = 10, l = 20, unit = "pt")
    )
}

save_pdf_jpg <- function(plot_obj, pdf_file, jpg_file, w = 13, h = 10, dpi = 900) {
  pdf(pdf_file, width = w, height = h)
  print(plot_obj)
  dev.off()

  jpeg(jpg_file, width = dpi * w, height = dpi * h, res = dpi)
  print(plot_obj)
  dev.off()
}

# Step 1: Read and auto-detect columns
res_df <- read_csv(deg_file, show_col_types = FALSE)

gene_col  <- pick_col(res_df, c("Gene_Symbol", "SYMBOL", "Symbol", "gene", "Gene", "hgnc_symbol"))
padj_col  <- pick_col(res_df, c("padj", "adj.P.Val", "FDR", "qvalue", "adj_pval", "p_adj"))
logfc_col <- pick_col(res_df, c("log2FoldChange", "logFC", "log2FC", "LFC"))

if (is.na(gene_col) || is.na(padj_col) || is.na(logfc_col)) {
  stop("Could not auto-detect required columns. Available: ", paste(colnames(res_df), collapse = ", "))
}

res_df2 <- res_df %>%
  mutate(
    gene_symbol    = .data[[gene_col]],
    padj           = as.numeric(.data[[padj_col]]),
    log2FoldChange = as.numeric(.data[[logfc_col]])
  )

# Step 2: Filter DEGs
deg_filtered <- res_df2 %>%
  filter(!is.na(padj), padj < padj_cutoff, abs(log2FoldChange) > logfc_cutoff) %>%
  filter(!is.na(gene_symbol), gene_symbol != "")

gene_symbols <- unique(deg_filtered$gene_symbol)
message(sprintf("Significant DEGs for enrichment: %d", length(gene_symbols)))

# Step 3: Convert SYMBOL to ENTREZ ID
entrez_ids <- bitr(
  gene_symbols,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
) %>%
  filter(!is.na(ENTREZID)) %>%
  distinct(ENTREZID, .keep_all = TRUE)

gene_list <- entrez_ids$ENTREZID
message(sprintf("Mapped to Entrez IDs: %d", length(gene_list)))

if (length(gene_list) < 5) {
  stop("Too few Entrez IDs mapped. Check gene symbols.")
}

# Step 4: KEGG enrichment
kegg_res <- enrichKEGG(
  gene = gene_list,
  organism = "hsa",
  pvalueCutoff = 0.05
)

# Remove COVID-19 terms if present
if (!is.null(kegg_res) && nrow(as.data.frame(kegg_res)) > 0) {
  kk <- as.data.frame(kegg_res)
  kk <- kk %>% filter(!grepl("COVID-19", Description, ignore.case = TRUE))
  kegg_res@result <- kk
}

# Step 5: GO enrichment (BP, CC, MF)
go_bp <- enrichGO(gene = gene_list, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                  ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

go_cc <- enrichGO(gene = gene_list, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                  ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

go_mf <- enrichGO(gene = gene_list, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                  ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

# Step 6: Save result CSVs
write.csv(as.data.frame(kegg_res), sprintf("results/KEGG_results_%s.csv", out_suffix), row.names = FALSE)
write.csv(as.data.frame(go_bp),    sprintf("results/GO_BP_results_%s.csv", out_suffix), row.names = FALSE)
write.csv(as.data.frame(go_cc),    sprintf("results/GO_CC_results_%s.csv", out_suffix), row.names = FALSE)
write.csv(as.data.frame(go_mf),    sprintf("results/GO_MF_results_%s.csv", out_suffix), row.names = FALSE)
write.csv(entrez_ids,              sprintf("results/SYMBOL_to_ENTREZID_%s.csv", out_suffix), row.names = FALSE)

# Step 7: Generate plots
if (!is.null(kegg_res) && nrow(as.data.frame(kegg_res)) > 0) {
  p1 <- dotplot(kegg_res, showCategory = show_kegg, label_format = 65) +
    ggtitle(sprintf("KEGG Pathway Enrichment - %s", out_suffix))
  p1 <- apply_readable_theme(p1)
  save_pdf_jpg(p1,
    sprintf("results/KEGG_dotplot_%s.pdf", out_suffix),
    sprintf("results/KEGG_dotplot_%s.jpg", out_suffix))

  p1b <- barplot(kegg_res, showCategory = show_kegg, x = "GeneRatio") +
    ggtitle(sprintf("KEGG Pathways (Barplot) - %s", out_suffix))
  p1b <- apply_readable_theme(p1b)
  save_pdf_jpg(p1b,
    sprintf("results/KEGG_barplot_%s.pdf", out_suffix),
    sprintf("results/KEGG_barplot_%s.jpg", out_suffix))
} else {
  warning("KEGG enrichment returned 0 terms. Skipping plots.")
}

if (!is.null(go_bp) && nrow(as.data.frame(go_bp)) > 0) {
  p2 <- dotplot(go_bp, showCategory = show_go, label_format = 65) +
    ggtitle(sprintf("GO Biological Process - %s", out_suffix))
  p2 <- apply_readable_theme(p2)
  save_pdf_jpg(p2,
    sprintf("results/GO_BP_dotplot_%s.pdf", out_suffix),
    sprintf("results/GO_BP_dotplot_%s.jpg", out_suffix))
}

if (!is.null(go_cc) && nrow(as.data.frame(go_cc)) > 0) {
  p3 <- dotplot(go_cc, showCategory = show_go, label_format = 65) +
    ggtitle(sprintf("GO Cellular Component - %s", out_suffix))
  p3 <- apply_readable_theme(p3)
  save_pdf_jpg(p3,
    sprintf("results/GO_CC_dotplot_%s.pdf", out_suffix),
    sprintf("results/GO_CC_dotplot_%s.jpg", out_suffix))
}

if (!is.null(go_mf) && nrow(as.data.frame(go_mf)) > 0) {
  p4 <- dotplot(go_mf, showCategory = show_go, label_format = 65) +
    ggtitle(sprintf("GO Molecular Function - %s", out_suffix))
  p4 <- apply_readable_theme(p4)
  save_pdf_jpg(p4,
    sprintf("results/GO_MF_dotplot_%s.pdf", out_suffix),
    sprintf("results/GO_MF_dotplot_%s.jpg", out_suffix))
}

message(sprintf("KEGG/GO enrichment complete (%s). Results saved in results/", out_suffix))
