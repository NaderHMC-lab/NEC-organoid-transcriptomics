#!/usr/bin/env python3
"""
Script: 05_deg_overlap.py
Project: Temporal Gene Expression Dynamics of NEC in Human Intestinal
         Organoid Model (Ellaithi, Gattu Linga et al.)

Purpose:
    Temporal overlap analysis of DEGs between 24 h and 72 h post-exposure
    to identify stage-specific and persistent transcriptional programs
    (Section 3.5 and Tables S5-S8 in the manuscript).

    Categorizes genes as:
      - Persistently upregulated at both 24 h and 72 h (Table S5)
      - Persistently downregulated at both 24 h and 72 h (Table S6)
      - Transiently upregulated at 24 h, then downregulated at 72 h (Table S7)
      - Delayed response: downregulated at 24 h, upregulated at 72 h (Table S8)

    Outputs CSVs with logFC and adjusted p-values from both timepoints,
    plus an optional consolidated Excel workbook.

Usage:
    python 05_deg_overlap.py \
        --deg24 results/LIMMA_DEG_significant_NEC_vs_Control_24h.csv \
        --deg72 results/LIMMA_DEG_significant_NEC_vs_Control_72h.csv \
        --outdir results/deg_overlaps \
        --xlsx
"""

import argparse
import pandas as pd
import numpy as np
from pathlib import Path


def read_any(path: Path) -> pd.DataFrame:
    """Read CSV, TSV, or Excel file."""
    try:
        return pd.read_csv(path, sep=None, engine="python")
    except Exception:
        try:
            return pd.read_csv(path, sep="\t")
        except Exception:
            return pd.read_excel(path)


def guess_col(cols, candidates, loose_keys=None):
    """Auto-detect column name from a list of candidates."""
    cset = {c.lower(): c for c in cols}
    for cand in candidates:
        if cand.lower() in cset:
            return cset[cand.lower()]
    if loose_keys:
        for c in cols:
            cl = c.lower()
            if all(k in cl for k in loose_keys):
                return c
    return None


def guess_gene_col(cols):
    return guess_col(cols,
        ["Gene_Symbol", "SYMBOL", "Symbol", "gene", "Gene", "hgnc_symbol", "gene_name"],
        ["symbol"])


def guess_logfc_col(cols):
    return guess_col(cols,
        ["logFC", "log2FC", "log2FoldChange", "LFC"],
        ["log", "fc"])


def guess_adjp_col(cols):
    return guess_col(cols,
        ["adj.P.Val", "adj_p_val", "padj", "FDR", "qvalue"],
        ["adj", "val"])


def clean_subset(df, gene_col, logfc_col, adjp_col):
    """Extract and clean relevant columns."""
    d = df[[gene_col, logfc_col, adjp_col]].copy()
    d = d.dropna(subset=[gene_col, logfc_col, adjp_col])
    d[gene_col] = d[gene_col].astype(str).str.strip()
    d = d[d[gene_col] != ""]
    d[logfc_col] = pd.to_numeric(d[logfc_col], errors="coerce")
    d[adjp_col] = pd.to_numeric(d[adjp_col], errors="coerce")
    d = d.dropna(subset=[logfc_col, adjp_col])
    return d


def dedup_keep_min_adjp(df, gene_col, logfc_col, adjp_col):
    """Keep one row per gene: smallest adj p-value, then largest |logFC|."""
    df = df.copy()
    df["_abslog"] = df[logfc_col].abs()
    df = df.sort_values([adjp_col, "_abslog"], ascending=[True, False])
    df = df.drop_duplicates(subset=[gene_col], keep="first")
    df = df.drop(columns=["_abslog"])
    return df


def main():
    ap = argparse.ArgumentParser(
        description="Find DEG overlaps between 24h and 72h timepoints.")
    ap.add_argument("--deg24", required=True, type=Path, help="24h DEG file")
    ap.add_argument("--deg72", required=True, type=Path, help="72h DEG file")
    ap.add_argument("--outdir", default="results/deg_overlaps", type=Path)
    ap.add_argument("--xlsx", action="store_true", help="Also write Excel workbook")
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    # Load data
    df24 = read_any(args.deg24)
    df72 = read_any(args.deg72)

    # Auto-detect columns
    gene24, lfc24, ap24 = guess_gene_col(df24.columns), guess_logfc_col(df24.columns), guess_adjp_col(df24.columns)
    gene72, lfc72, ap72 = guess_gene_col(df72.columns), guess_logfc_col(df72.columns), guess_adjp_col(df72.columns)

    if not all([gene24, lfc24, ap24, gene72, lfc72, ap72]):
        raise SystemExit(
            f"Could not detect required columns.\n"
            f"24h: gene='{gene24}' logFC='{lfc24}' adjP='{ap24}'\n"
            f"72h: gene='{gene72}' logFC='{lfc72}' adjP='{ap72}'"
        )

    # Clean, rename, deduplicate
    d24 = clean_subset(df24, gene24, lfc24, ap24).rename(
        columns={gene24: "Gene", lfc24: "logFC_24h", ap24: "adjP_24h"})
    d72 = clean_subset(df72, gene72, lfc72, ap72).rename(
        columns={gene72: "Gene", lfc72: "logFC_72h", ap72: "adjP_72h"})

    d24 = dedup_keep_min_adjp(d24, "Gene", "logFC_24h", "adjP_24h")
    d72 = dedup_keep_min_adjp(d72, "Gene", "logFC_72h", "adjP_72h")

    # Assign direction
    d24["dir_24h"] = np.where(d24["logFC_24h"] > 0, "up", np.where(d24["logFC_24h"] < 0, "down", "zero"))
    d72["dir_72h"] = np.where(d72["logFC_72h"] > 0, "up", np.where(d72["logFC_72h"] < 0, "down", "zero"))

    d24 = d24[d24["dir_24h"] != "zero"]
    d72 = d72[d72["dir_72h"] != "zero"]

    # Split by direction
    d24u = d24[d24["dir_24h"] == "up"].set_index("Gene")
    d24d = d24[d24["dir_24h"] == "down"].set_index("Gene")
    d72u = d72[d72["dir_72h"] == "up"].set_index("Gene")
    d72d = d72[d72["dir_72h"] == "down"].set_index("Gene")

    # Find intersections
    common_up = sorted(set(d24u.index) & set(d72u.index))
    common_down = sorted(set(d24d.index) & set(d72d.index))
    up24_down72 = sorted(set(d24u.index) & set(d72d.index))
    down24_up72 = sorted(set(d24d.index) & set(d72u.index))

    # Build result tables
    common_up_df = (d24u[["logFC_24h", "adjP_24h"]]
                    .join(d72u[["logFC_72h", "adjP_72h"]], how="inner")
                    ).loc[common_up].reset_index()

    common_down_df = (d24d[["logFC_24h", "adjP_24h"]]
                      .join(d72d[["logFC_72h", "adjP_72h"]], how="inner")
                      ).loc[common_down].reset_index()

    discord_ud = (d24u[["logFC_24h", "adjP_24h"]]
                  .join(d72d[["logFC_72h", "adjP_72h"]], how="inner")
                  ).loc[up24_down72].reset_index()

    discord_du = (d24d[["logFC_24h", "adjP_24h"]]
                  .join(d72u[["logFC_72h", "adjP_72h"]], how="inner")
                  ).loc[down24_up72].reset_index()

    # Save CSVs
    common_up_df.to_csv(args.outdir / "common_up_24h_72h.csv", index=False)
    common_down_df.to_csv(args.outdir / "common_down_24h_72h.csv", index=False)
    discord_ud.to_csv(args.outdir / "discordant_up24_down72.csv", index=False)
    discord_du.to_csv(args.outdir / "discordant_down24_up72.csv", index=False)

    summary = pd.DataFrame({
        "Category": ["24h up", "24h down", "72h up", "72h down",
                     "Common up", "Common down", "Up24-Down72", "Down24-Up72"],
        "Count": [len(d24u), len(d24d), len(d72u), len(d72d),
                  len(common_up), len(common_down),
                  len(up24_down72), len(down24_up72)]
    })
    summary.to_csv(args.outdir / "summary_counts.tsv", sep="\t", index=False)

    # Optional Excel
    if args.xlsx:
        xlsx_path = args.outdir / "deg_overlaps.xlsx"
        with pd.ExcelWriter(xlsx_path, engine="xlsxwriter") as xw:
            summary.to_excel(xw, index=False, sheet_name="Summary")
            common_up_df.to_excel(xw, index=False, sheet_name="Common_Up")
            common_down_df.to_excel(xw, index=False, sheet_name="Common_Down")
            discord_ud.to_excel(xw, index=False, sheet_name="Up24_Down72")
            discord_du.to_excel(xw, index=False, sheet_name="Down24_Up72")

    print("=== DEG Overlap Summary ===")
    print(summary.to_string(index=False))
    print(f"\nResults saved to: {args.outdir}/")


if __name__ == "__main__":
    main()
