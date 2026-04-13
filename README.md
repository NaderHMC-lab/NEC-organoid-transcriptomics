# Temporal Gene Expression Dynamics of Necrotizing Enterocolitis (NEC) in Human Intestinal Organoid Model

## Overview

This repository contains the analysis scripts for the study:

> **Temporal Gene Expression Dynamics of Necrotizing Enterocolitis (NEC) in Human Intestinal Organoid Model**
>
> Mona Ellaithi, BalaSubramani Gattu Linga, Waleed Aamer, Faisal E. Ibrahim, Jameela Roshanuddin, Rand Hamdan, Nader Al-Dewik

We performed time-resolved reanalysis of RNA-seq data from neonatal human intestinal organoids exposed to NEC-associated bacteria ([GEO: GSE226086](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE226086), BioProject: [PRJNA938490](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA938490)). Differential expression at 24 h and 72 h was identified using LIMMA (|log2FC| >= 1, adj P < 0.05), followed by KEGG, Gene Ontology, and temporal overlap analyses to define stage-specific and persistent transcriptional programs associated with NEC progression.

### Study Design

The dataset was generated using a neonatal human intestinal organoid-based microfluidic platform. Neonatal small intestinal epithelial organoids were cultured within microfluidic chips and exposed to bacterial isolates derived from preterm infants diagnosed with NEC. RNA was collected at two timepoints:

| Group | Timepoint | Samples (n) |
|-------|-----------|-------------|
| Control (uninfected) | 24 h | 3 |
| NEC-exposed | 24 h | 3 |
| Control (uninfected) | 72 h | 3 |
| NEC-exposed | 72 h | 3 |

### Key Findings

- **24 h:** Early inflammatory and epithelial stress responses with upregulation of immune (CSF2, CSF3, VCAM1), translational (EEF1G), and metabolic genes (SI, APOA4), and suppression of cell-cycle (MCM4, BIRC5, BRCA2) and epithelial transport (SLC2A2, SLC4A8) pathways.
- **72 h:** Intensified transcriptional dysregulation with mitochondrial hyperactivation (MT-ND5, MT-CO1, MT-CYB), metabolic reprogramming (APOA1, APOC3, FTL), and broad repression of chromatin regulators (SPEN, ASH1L, TET3), DNA damage-response genes (PRKDC, SETD2), and cell-cycle machinery.
- **Temporal overlap:** 271 genes persistently upregulated and 2,035 genes persistently downregulated at both timepoints, with systematic amplification of effect sizes from 24 h to 72 h, indicating consolidation of pathogenic transcriptional programs.

## Analysis Pipeline

The pipeline consists of six scripts executed sequentially. Scripts 01-04 (R) are run separately for each timepoint by setting a `timepoint` configuration variable. Scripts 05-06 (Python) integrate results across both timepoints.

```
01_limma_deg_analysis.R        Differential expression (LIMMA)
        |
        v
02_volcano_plot.R              Volcano plot visualization
03_heatmap_top50_degs.R        Top 50 DEG heatmap
04_kegg_go_enrichment.R        KEGG + GO enrichment
        |
        v
05_deg_overlap.py              Temporal overlap analysis (24h vs 72h)
06_export_top_common_degs.py   Export top common DEGs
```

### Script Descriptions

| Script | Description |
|--------|-------------|
| `01_limma_deg_analysis.R` | Reads normalized gene expression data, removes lowly expressed genes (mean < 1), applies log2(x + 1) transformation, fits LIMMA linear models with empirical Bayes moderation, and extracts significant DEGs (\|log2FC\| >= 1, BH-adjusted P < 0.05). |
| `02_volcano_plot.R` | Generates publication-quality volcano plots showing the global distribution of DEGs with significance thresholds. Top 5 upregulated and top 5 downregulated genes are labeled. Output at 600 DPI. |
| `03_heatmap_top50_degs.R` | Selects the top 25 upregulated and 25 downregulated DEGs by log2FC, extracts their expression profiles from the normalized matrix, applies Z-score normalization, and generates annotated heatmaps with hierarchical clustering (PNG + PDF). |
| `04_kegg_go_enrichment.R` | Converts gene symbols to Entrez IDs (human, org.Hs.eg.db) and performs over-representation analysis for KEGG pathways and GO terms (BP, CC, MF). Generates dot plots and bar plots (PDF + 900 DPI JPG). COVID-19-related KEGG terms are excluded. |
| `05_deg_overlap.py` | Compares DEG lists from 24 h and 72 h to identify commonly upregulated, commonly downregulated, and discordant (direction-reversed) genes across timepoints. Outputs CSVs and optional Excel workbook with summary counts. |
| `06_export_top_common_degs.py` | Extracts the top N (default 10) common upregulated and downregulated genes between 24 h and 72 h, ranked by mean log2FC, with p-values and adjusted p-values from both timepoints. Outputs a single Excel workbook. |

## Directory Structure

```
NEC_scripts/
├── data/                              # Input data (not tracked in git)
│   ├── GSE226086_normalized_gene_expression_24hr.xlsx
│   └── GSE226086_normalized_gene_expression_72hr.xlsx
├── results/                           # All output files (not tracked in git)
│   ├── LIMMA_DEG_significant_NEC_vs_Control_24h.csv
│   ├── LIMMA_DEG_significant_NEC_vs_Control_72h.csv
│   ├── volcano_plot_NEC_vs_Control_*.png
│   ├── NEC_vs_Control_*_heatmap_top50.{png,pdf}
│   ├── KEGG_*.{csv,pdf,jpg}
│   ├── GO_*_*.{csv,pdf,jpg}
│   ├── deg_overlaps/
│   └── top10_common_degs.xlsx
├── 01_limma_deg_analysis.R
├── 02_volcano_plot.R
├── 03_heatmap_top50_degs.R
├── 04_kegg_go_enrichment.R
├── 05_deg_overlap.py
├── 06_export_top_common_degs.py
├── run_pipeline.sh                    # Run full pipeline (both timepoints)
├── install_R_packages.R               # One-click R dependency installer
├── requirements.txt                   # Python dependencies
├── CITATION.cff                       # Machine-readable citation metadata
├── LICENSE                            # MIT License
├── .gitignore
├── .github/
│   └── ISSUE_TEMPLATE/
│       ├── bug_report.md
│       └── question.md
└── README.md
```

## Quick Start

```bash
# Clone the repository
git clone https://github.com/glbala87/NEC-organoid-transcriptomics.git
cd NEC-organoid-transcriptomics

# Install dependencies
Rscript install_R_packages.R
pip install -r requirements.txt

# Place input data in data/ directory (download from GEO: GSE226086)

# Run the full pipeline
chmod +x run_pipeline.sh
./run_pipeline.sh
```

## Requirements

### Software

- **R** >= 4.5.0
- **Python** >= 3.8

### R Packages

```bash
Rscript install_R_packages.R
```

Or install manually:
```r
# CRAN
install.packages(c("readxl", "dplyr", "ggplot2", "ggrepel", "pheatmap",
                   "RColorBrewer", "tibble", "readr", "stringr"))

# Bioconductor
BiocManager::install(c("limma", "clusterProfiler", "org.Hs.eg.db", "enrichplot"))
```

### Python Packages

```bash
pip install -r requirements.txt
```

## Usage

### Option A: Run Full Pipeline (Recommended)

```bash
./run_pipeline.sh
```

This runs all scripts for both timepoints automatically and reports results.

### Option B: Run Scripts Individually

**Step 1: Prepare Input Data**

Download the normalized gene expression file from [GEO: GSE226086](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE226086) and place the Excel files in the `data/` directory:
- `GSE226086_normalized_gene_expression_24hr.xlsx`
- `GSE226086_normalized_gene_expression_72hr.xlsx`

**Step 2: Run DEG Analysis (R)**

Each R script contains a `timepoint` variable at the top. Set it to `"24h"` or `"72h"` and source the script. Run each timepoint separately:

```r
# Set working directory to NEC-organoid-transcriptomics/
setwd("path/to/NEC-organoid-transcriptomics")

# --- For 24 h analysis ---
# Edit timepoint <- "24h" in each script, then:
source("01_limma_deg_analysis.R")
source("02_volcano_plot.R")
source("03_heatmap_top50_degs.R")
source("04_kegg_go_enrichment.R")

# --- For 72 h analysis ---
# Edit timepoint <- "72h" in each script, then:
source("01_limma_deg_analysis.R")
source("02_volcano_plot.R")
source("03_heatmap_top50_degs.R")
source("04_kegg_go_enrichment.R")
```

**Step 3: Temporal Overlap Analysis (Python)**

After running both timepoints through scripts 01-04:

```bash
# Identify common and discordant DEGs between 24h and 72h
python 05_deg_overlap.py \
  --deg24 results/LIMMA_DEG_significant_NEC_vs_Control_24h.csv \
  --deg72 results/LIMMA_DEG_significant_NEC_vs_Control_72h.csv \
  --outdir results/deg_overlaps \
  --xlsx

# Export top 10 common DEGs ranked by mean logFC
python 06_export_top_common_degs.py \
  --deg24 results/LIMMA_DEG_significant_NEC_vs_Control_24h.csv \
  --deg72 results/LIMMA_DEG_significant_NEC_vs_Control_72h.csv \
  --out results/top10_common_degs.xlsx
```

## Analysis Parameters

| Parameter | Value | Script |
|-----------|-------|--------|
| Low expression filter | Mean normalized expression >= 1 | 01 |
| Log2 transformation | log2(x + 1) | 01 |
| Differential expression method | LIMMA (lmFit + eBayes) | 01 |
| DEG significance thresholds | \|log2FC\| >= 1, BH-adjusted P < 0.05 | 01 |
| P-value adjustment | Benjamini-Hochberg (BH) | 01, 04 |
| Top labeled genes (volcano) | 5 upregulated + 5 downregulated | 02 |
| Top DEGs for heatmap | 25 upregulated + 25 downregulated | 03 |
| Heatmap normalization | Row-wise Z-score | 03 |
| Enrichment p-value cutoff | 0.05 | 04 |
| Organism database | Human (hsa / org.Hs.eg.db) | 04 |
| Enrichment categories shown | Top 20 per category | 04 |
| Plot resolution | 600 DPI (volcano, heatmap), 900 DPI (enrichment) | 02-04 |

## Corresponding Manuscript Figures

| Figure | Generated by | Description |
|--------|-------------|-------------|
| Figure 1a | `02_volcano_plot.R` (24h) | Volcano plot of DEGs at 24 h |
| Figure 1b | `03_heatmap_top50_degs.R` (24h) | Heatmap of top 50 DEGs at 24 h |
| Figure 2a | `02_volcano_plot.R` (72h) | Volcano plot of DEGs at 72 h |
| Figure 2b | `03_heatmap_top50_degs.R` (72h) | Heatmap of top 50 DEGs at 72 h |
| Figure 3 | `04_kegg_go_enrichment.R` (24h) | KEGG pathway enrichment at 24 h |
| Figure 4 | `04_kegg_go_enrichment.R` (24h) | GO enrichment (BP, CC, MF) at 24 h |
| Figure 5 | `04_kegg_go_enrichment.R` (72h) | KEGG pathway enrichment at 72 h |
| Figure 6 | `04_kegg_go_enrichment.R` (72h) | GO enrichment (BP, CC, MF) at 72 h |
| Tables S1-S4 | `01_limma_deg_analysis.R` | Significant DEGs at 24 h and 72 h |
| Tables S5-S8 | `05_deg_overlap.py` | Temporal overlap: common and discordant DEGs |

## Data Availability

The dataset analyzed in this study was obtained from the Gene Expression Omnibus (NCBI) under accession number [GSE226086](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE226086) (BioProject: [PRJNA938490](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA938490)).

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

## Contact

For questions or data requests, contact the corresponding author:

**Nader Al-Dewik**
- Translational and Precision Medicine Research, Hamad Medical Corporation, Qatar
- Email: naldewik@hamad.qa / nader.al-dewik@kingston.ac.uk
