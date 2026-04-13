#!/usr/bin/env bash
# ==============================================================================
# Script: run_pipeline.sh
# Purpose: Run the complete NEC transcriptomics DEG analysis pipeline for
#          both timepoints (24h and 72h).
#
# Prerequisites:
#   - R >= 4.5.0 with all packages installed (run: Rscript install_R_packages.R)
#   - Python >= 3.8 with packages installed (run: pip install -r requirements.txt)
#   - Input data files in data/ directory
#
# Usage:
#   chmod +x run_pipeline.sh
#   ./run_pipeline.sh
# ==============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Create output directories
mkdir -p results/deg_overlaps

# ------------------------------------------------------------------------------
# Helper function
# ------------------------------------------------------------------------------
run_step() {
    local step="$1"
    local desc="$2"
    echo ""
    echo "================================================================"
    echo "  $desc"
    echo "================================================================"
    eval "$step"
    echo "  Done."
}

# ------------------------------------------------------------------------------
# Check input data
# ------------------------------------------------------------------------------
echo "Checking input data..."
for f in data/GSE226086_normalized_gene_expression_24hr.xlsx \
         data/GSE226086_normalized_gene_expression_72hr.xlsx; do
    if [ ! -f "$f" ]; then
        echo "ERROR: Missing input file: $f"
        echo "Download from GEO (GSE226086) and place in data/ directory."
        exit 1
    fi
done
echo "Input data found."

# ==============================================================================
# 24-HOUR ANALYSIS
# ==============================================================================
echo ""
echo "╔══════════════════════════════════════════════════════════════════╗"
echo "║                    24-HOUR TIMEPOINT                           ║"
echo "╚══════════════════════════════════════════════════════════════════╝"

# Set timepoint to 24h in R scripts using temp modified copies
export NEC_TIMEPOINT="24h"

run_step \
    "Rscript -e 'timepoint <- \"24h\"; source(\"01_limma_deg_analysis.R\", local=TRUE)'" \
    "Step 1/4: LIMMA DEG analysis (24h)"

run_step \
    "Rscript -e 'timepoint <- \"24h\"; source(\"02_volcano_plot.R\", local=TRUE)'" \
    "Step 2/4: Volcano plot (24h)"

run_step \
    "Rscript -e 'timepoint <- \"24h\"; source(\"03_heatmap_top50_degs.R\", local=TRUE)'" \
    "Step 3/4: Heatmap of top 50 DEGs (24h)"

run_step \
    "Rscript -e 'timepoint <- \"24h\"; source(\"04_kegg_go_enrichment.R\", local=TRUE)'" \
    "Step 4/4: KEGG and GO enrichment (24h)"

# ==============================================================================
# 72-HOUR ANALYSIS
# ==============================================================================
echo ""
echo "╔══════════════════════════════════════════════════════════════════╗"
echo "║                    72-HOUR TIMEPOINT                           ║"
echo "╚══════════════════════════════════════════════════════════════════╝"

export NEC_TIMEPOINT="72h"

run_step \
    "Rscript -e 'timepoint <- \"72h\"; source(\"01_limma_deg_analysis.R\", local=TRUE)'" \
    "Step 1/4: LIMMA DEG analysis (72h)"

run_step \
    "Rscript -e 'timepoint <- \"72h\"; source(\"02_volcano_plot.R\", local=TRUE)'" \
    "Step 2/4: Volcano plot (72h)"

run_step \
    "Rscript -e 'timepoint <- \"72h\"; source(\"03_heatmap_top50_degs.R\", local=TRUE)'" \
    "Step 3/4: Heatmap of top 50 DEGs (72h)"

run_step \
    "Rscript -e 'timepoint <- \"72h\"; source(\"04_kegg_go_enrichment.R\", local=TRUE)'" \
    "Step 4/4: KEGG and GO enrichment (72h)"

# ==============================================================================
# CROSS-TIMEPOINT ANALYSIS
# ==============================================================================
echo ""
echo "╔══════════════════════════════════════════════════════════════════╗"
echo "║               CROSS-TIMEPOINT ANALYSIS                        ║"
echo "╚══════════════════════════════════════════════════════════════════╝"

run_step \
    "python3 05_deg_overlap.py \
        --deg24 results/LIMMA_DEG_significant_NEC_vs_Control_24h.csv \
        --deg72 results/LIMMA_DEG_significant_NEC_vs_Control_72h.csv \
        --outdir results/deg_overlaps \
        --xlsx" \
    "Step 1/2: DEG temporal overlap analysis"

run_step \
    "python3 06_export_top_common_degs.py \
        --deg24 results/LIMMA_DEG_significant_NEC_vs_Control_24h.csv \
        --deg72 results/LIMMA_DEG_significant_NEC_vs_Control_72h.csv \
        --out results/top10_common_degs.xlsx" \
    "Step 2/2: Export top common DEGs"

# ==============================================================================
# SUMMARY
# ==============================================================================
echo ""
echo "╔══════════════════════════════════════════════════════════════════╗"
echo "║                  PIPELINE COMPLETE                             ║"
echo "╚══════════════════════════════════════════════════════════════════╝"
echo ""
echo "All results saved in: results/"
echo ""
echo "Output files:"
find results/ -type f | sort | sed 's/^/  /'
