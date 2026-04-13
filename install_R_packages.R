# ==============================================================================
# Script: install_R_packages.R
# Purpose: Install all R package dependencies for the NEC transcriptomics
#          analysis pipeline.
#
# Usage:
#   Rscript install_R_packages.R
#   # or source("install_R_packages.R") from within R/RStudio
# ==============================================================================

cat("Installing CRAN packages...\n")

cran_packages <- c(
  "readxl",        # Read Excel files
  "dplyr",         # Data manipulation
  "ggplot2",       # Plotting
  "ggrepel",       # Non-overlapping text labels
  "pheatmap",      # Heatmap generation
  "RColorBrewer",  # Color palettes
  "tibble",        # Tibble data structures
  "readr",         # Read CSV/TSV files
  "stringr"        # String manipulation
)

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  Installing %s...\n", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  } else {
    cat(sprintf("  %s already installed\n", pkg))
  }
}

cat("\nInstalling Bioconductor packages...\n")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

bioc_packages <- c(
  "limma",            # Linear models for microarray/RNA-seq
  "clusterProfiler",  # Enrichment analysis
  "org.Hs.eg.db",     # Human gene annotation
  "enrichplot"        # Enrichment visualization
)

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  Installing %s...\n", pkg))
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  } else {
    cat(sprintf("  %s already installed\n", pkg))
  }
}

cat("\nAll dependencies installed successfully.\n")
cat("Run sessionInfo() to verify package versions.\n")
