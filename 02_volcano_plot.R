# ==============================================================================
# Script: 02_volcano_plot.R
# Project: Temporal Gene Expression Dynamics of NEC in Human Intestinal
#          Organoid Model (Ellaithi, Gattu Linga et al.)
#
# Purpose:
#   Generate publication-quality volcano plots of LIMMA differential expression
#   results (Figures 1a and 2a in the manuscript).
#
# Method:
#   - Plots log2 fold change (x-axis) vs -log10 adjusted P-value (y-axis)
#   - Applies significance thresholds: |log2FC| >= 1, adj.P.Val < 0.05
#   - Color-codes genes as Upregulated (red), Downregulated (blue), or
#     Not significant (gray)
#   - Labels the top 5 upregulated and top 5 downregulated genes by adj.P.Val
#
# Input:  results/LIMMA_DEG_significant_NEC_vs_Control_{24h,72h}.csv
# Output: results/volcano_plot_NEC_vs_Control_{24h,72h}.png (600 DPI)
#
# Usage:
#   Set `timepoint` below to "24h" or "72h", then source this script.
# ==============================================================================

# ---- Configuration -----------------------------------------------------------
timepoint <- "24h"  # Change to "72h" for 72-hour analysis
# ------------------------------------------------------------------------------

library(ggplot2)
library(ggrepel)
library(dplyr)

# Define timepoint-specific parameters
if (timepoint == "24h") {
  input_file  <- "results/LIMMA_DEG_significant_NEC_vs_Control_24h.csv"
  output_file <- "results/volcano_plot_NEC_vs_Control_24h.png"
  plot_title  <- "Volcano Plot: NEC vs Control - 24h (LIMMA)"
} else if (timepoint == "72h") {
  input_file  <- "results/LIMMA_DEG_significant_NEC_vs_Control_72h.csv"
  output_file <- "results/volcano_plot_NEC_vs_Control_72h.png"
  plot_title  <- "Volcano Plot: NEC vs Control - 72h (LIMMA)"
} else {
  stop("timepoint must be '24h' or '72h'")
}

# Step 1: Load DEG results
deg_df <- read.csv(input_file)

# Step 2: Annotate significance
deg_df <- deg_df %>%
  mutate(
    neg_log10_adjP = -log10(adj.P.Val),
    significance = case_when(
      adj.P.Val < 0.05 & logFC >= 1  ~ "Upregulated",
      adj.P.Val < 0.05 & logFC <= -1 ~ "Downregulated",
      TRUE ~ "Not significant"
    )
  )

# Step 3: Select top genes for labeling
top_up   <- deg_df %>% filter(significance == "Upregulated") %>% arrange(adj.P.Val) %>% slice(1:5)
top_down <- deg_df %>% filter(significance == "Downregulated") %>% arrange(adj.P.Val) %>% slice(1:5)
top_genes <- bind_rows(top_up, top_down)

deg_df$label <- ifelse(deg_df$Gene_Symbol %in% top_genes$Gene_Symbol,
                       deg_df$Gene_Symbol, NA)

# Step 4: Generate volcano plot
volcano_plot <- ggplot(deg_df, aes(x = logFC, y = neg_log10_adjP)) +
  geom_point(aes(color = significance), size = 2, alpha = 0.7, stroke = 0.1) +
  scale_color_manual(values = c(
    "Upregulated"     = "#D62728",
    "Downregulated"   = "#1F77B4",
    "Not significant" = "gray70"
  )) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_text_repel(
    aes(label = label),
    size = 4,
    box.padding = 0.3,
    point.padding = 0.2,
    max.overlaps = 100,
    segment.size = 0.3,
    segment.alpha = 0.7,
    min.segment.length = 0
  ) +
  labs(
    title = plot_title,
    x = expression(log[2]~Fold~Change),
    y = expression(-log[10]~Adjusted~P~Value),
    color = "Regulation"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

# Step 5: Save plot
ggsave(output_file, plot = volcano_plot, width = 8, height = 6, dpi = 600)

cat(sprintf("Volcano plot saved: %s\n", output_file))
