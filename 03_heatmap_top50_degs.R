# ==============================================================================
# Script: 03_heatmap_top50_degs.R
# Project: Temporal Gene Expression Dynamics of NEC in Human Intestinal
#          Organoid Model (Ellaithi, Gattu Linga et al.)
#
# Purpose:
#   Generate hierarchically clustered heatmaps of the top 50 differentially
#   expressed genes (25 upregulated + 25 downregulated) between NEC-exposed
#   and control intestinal organoids (Figures 1b and 2b in the manuscript).
#
# Method:
#   - Selects top 25 upregulated and top 25 downregulated DEGs by log2FC
#   - Extracts their expression profiles from the normalized expression matrix
#   - Applies row-wise Z-score normalization to emphasize relative differences
#   - Generates annotated heatmaps with hierarchical row clustering
#   - Columns ordered as Control (n=3) then NEC (n=3), not clustered
#
# Input:
#   - results/LIMMA_DEG_significant_NEC_vs_Control_{24h,72h}.csv
#   - data/GSE226086_normalized_gene_expression_{24hr,72hr}.xlsx
# Output:
#   - results/NEC_vs_Control_{24h,72h}_top25_upregulated.csv
#   - results/NEC_vs_Control_{24h,72h}_top25_downregulated.csv
#   - results/NEC_vs_Control_{24h,72h}_top50_DEGs.csv
#   - results/NEC_vs_Control_{24h,72h}_heatmap_top50.{png,pdf}
#
# Usage:
#   Set `timepoint` below to "24h" or "72h", then source this script.
# ==============================================================================

# ---- Configuration -----------------------------------------------------------
timepoint <- "24h"  # Change to "72h" for 72-hour analysis
# ------------------------------------------------------------------------------

library(readxl)
library(pheatmap)
library(dplyr)
library(RColorBrewer)
library(tibble)

# Define timepoint-specific parameters
if (timepoint == "24h") {
  deg_file    <- "results/LIMMA_DEG_significant_NEC_vs_Control_24h.csv"
  expr_file   <- "data/GSE226086_normalized_gene_expression_24hr.xlsx"
  sample_order <- c("24_hr_control_1", "24_hr_control_2", "24_hr_control_3",
                    "24_hr_nec_1", "24_hr_nec_2", "24_hr_nec_3")
  out_prefix  <- "results/NEC_vs_Control_24h"
  plot_title  <- "Top 50 DEGs (25 Up, 25 Down) in NEC vs Control (24h)"
} else if (timepoint == "72h") {
  deg_file    <- "results/LIMMA_DEG_significant_NEC_vs_Control_72h.csv"
  expr_file   <- "data/GSE226086_normalized_gene_expression_72hr.xlsx"
  sample_order <- c("72_hr_control_1", "72_hr_control_2", "72_hr_control_3",
                    "72_hr_nec_1", "72_hr_nec_2", "72_hr_nec_3")
  out_prefix  <- "results/NEC_vs_Control_72h"
  plot_title  <- "Top 50 DEGs (25 Up, 25 Down) in NEC vs Control (72h)"
} else {
  stop("timepoint must be '24h' or '72h'")
}

# Step 1: Load DEG data and select top 25 up/down
deg_df <- read.csv(deg_file)

top_up   <- deg_df %>% arrange(desc(logFC)) %>% slice(1:25)
top_down <- deg_df %>% arrange(logFC) %>% slice(1:25)

write.csv(top_up,   paste0(out_prefix, "_top25_upregulated.csv"), row.names = FALSE)
write.csv(top_down, paste0(out_prefix, "_top25_downregulated.csv"), row.names = FALSE)

top50_combined <- bind_rows(top_up, top_down)
write.csv(top50_combined, paste0(out_prefix, "_top50_DEGs.csv"), row.names = FALSE)

# Step 2: Load and prepare expression matrix
expr_df <- read_excel(expr_file)
expr_df_clean <- expr_df %>%
  select(-gene_biotype) %>%
  distinct(Gene_Symbol, .keep_all = TRUE)

expr_mat <- expr_df_clean %>%
  column_to_rownames("Gene_Symbol")

# Step 3: Filter and Z-score normalize
top50_genes <- top50_combined$Gene_Symbol
expr_top50 <- expr_mat[rownames(expr_mat) %in% top50_genes, ]
expr_top50 <- expr_top50[match(top50_genes, rownames(expr_top50)), ]
expr_top50 <- expr_top50[, sample_order]

expr_scaled <- t(scale(t(expr_top50)))

# Step 4: Sample annotations
annotation_col <- data.frame(
  Group = factor(c(rep("Control", 3), rep("NEC", 3)),
                 levels = c("Control", "NEC"))
)
rownames(annotation_col) <- sample_order

# Step 5: Save heatmap (PNG)
png(paste0(out_prefix, "_heatmap_top50.png"), width = 2000, height = 2500, res = 300)
pheatmap(
  expr_scaled,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  show_rownames = TRUE,
  fontsize_row = 7,
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
  main = plot_title,
  border_color = NA
)
dev.off()

# Step 6: Save heatmap (PDF)
pdf(paste0(out_prefix, "_heatmap_top50.pdf"), width = 10, height = 12)
pheatmap(
  expr_scaled,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  show_rownames = TRUE,
  fontsize_row = 7,
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
  main = plot_title,
  border_color = NA
)
dev.off()

cat(sprintf("Heatmap saved: %s_heatmap_top50.{png,pdf}\n", out_prefix))
