# ==============================================================================
# Script: 01_limma_deg_analysis.R
# Project: Temporal Gene Expression Dynamics of NEC in Human Intestinal
#          Organoid Model (Ellaithi, Gattu Linga et al.)
# Dataset: GSE226086 (BioProject: PRJNA938490)
#
# Purpose:
#   Differential gene expression analysis comparing NEC-exposed vs control
#   neonatal human intestinal organoids at a specified timepoint (24h or 72h)
#   using the LIMMA framework.
#
# Method:
#   1. Load processed normalized gene expression matrix from GEO
#   2. Remove duplicate gene symbols (retain first occurrence)
#   3. Filter lowly expressed genes (mean normalized expression < 1)
#   4. Apply log2(x + 1) transformation (data is not pre-log-transformed)
#   5. Fit linear models (lmFit) with empirical Bayes moderation (eBayes)
#   6. Extract DEGs with |log2FC| >= 1 and BH-adjusted P < 0.05
#
# Input:  data/GSE226086_normalized_gene_expression_{24hr,72hr}.xlsx
# Output: results/LIMMA_DEG_significant_NEC_vs_Control_{24h,72h}.csv
#
# Usage:
#   Set `timepoint` below to "24h" or "72h", then source this script.
#
# Reference:
#   LIMMA: Ritchie et al., Nucleic Acids Res (2015)
# ==============================================================================

# ---- Configuration -----------------------------------------------------------
timepoint <- "24h"  # Change to "72h" for 72-hour analysis
# ------------------------------------------------------------------------------

library(readxl)
library(limma)
library(dplyr)

# Define timepoint-specific parameters
if (timepoint == "24h") {
  input_file     <- "data/GSE226086_normalized_gene_expression_24hr.xlsx"
  case_samples   <- c("24_hr_nec_1", "24_hr_nec_2", "24_hr_nec_3")
  control_samples <- c("24_hr_control_1", "24_hr_control_2", "24_hr_control_3")
  output_file    <- "results/LIMMA_DEG_significant_NEC_vs_Control_24h.csv"
} else if (timepoint == "72h") {
  input_file     <- "data/GSE226086_normalized_gene_expression_72hr.xlsx"
  case_samples   <- c("72_hr_nec_1", "72_hr_nec_2", "72_hr_nec_3")
  control_samples <- c("72_hr_control_1", "72_hr_control_2", "72_hr_control_3")
  output_file    <- "results/LIMMA_DEG_significant_NEC_vs_Control_72h.csv"
} else {
  stop("timepoint must be '24h' or '72h'")
}

# Step 1: Load and clean data
data <- read_excel(input_file)
data <- data %>% distinct(Gene_Symbol, .keep_all = TRUE)

gene_names <- data$Gene_Symbol
expression_data <- data %>% select(-Gene_Symbol, -gene_biotype)

# Step 2: Validate sample names
all_samples <- c(case_samples, control_samples)
stopifnot(all(all_samples %in% colnames(expression_data)))
expression_data <- expression_data[, all_samples]

# Step 3: Design matrix
group <- factor(c(rep("NEC", length(case_samples)),
                  rep("Control", length(control_samples))))
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# Step 4: Expression matrix
expression_matrix <- as.matrix(expression_data)
rownames(expression_matrix) <- gene_names

# Step 5: Filter lowly expressed genes (mean expression < 1)
row_means <- rowMeans(expression_matrix)
keep_genes <- row_means >= 1
expression_matrix <- expression_matrix[keep_genes, ]
gene_names <- gene_names[keep_genes]

# Step 6: Log2-transform
expression_matrix <- log2(expression_matrix + 1)

# Step 7: LIMMA linear model fitting
fit <- lmFit(expression_matrix, design)

contrast_matrix <- makeContrasts(NEC_vs_Control = NEC - Control, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Step 8: Extract and filter results
res <- topTable(fit2, coef = "NEC_vs_Control", adjust.method = "BH", number = Inf)
res$Gene_Symbol <- rownames(res)
res <- res %>% select(Gene_Symbol, logFC, P.Value, adj.P.Val, everything())

sig_DEGs <- res %>%
  filter(abs(logFC) > 1, adj.P.Val < 0.05)

# Step 9: Save results
write.csv(sig_DEGs, output_file, row.names = FALSE)

cat(sprintf("DEG analysis complete (%s): %d significant DEGs saved to %s\n",
            timepoint, nrow(sig_DEGs), output_file))
