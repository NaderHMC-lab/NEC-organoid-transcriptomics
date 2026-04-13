#!/usr/bin/env python3
"""
Script: 06_export_top_common_degs.py
Project: Temporal Gene Expression Dynamics of NEC in Human Intestinal
         Organoid Model (Ellaithi, Gattu Linga et al.)

Purpose:
    Export the top N persistently upregulated and persistently downregulated
    genes shared between 24 h and 72 h into a single Excel file for
    supplementary reporting and downstream interpretation.

    Genes are ranked by mean log2FC across both timepoints. The output
    includes logFC, P-value, and adjusted P-value from each timepoint.

Usage:
    python 06_export_top_common_degs.py \
        --deg24 results/LIMMA_DEG_significant_NEC_vs_Control_24h.csv \
        --deg72 results/LIMMA_DEG_significant_NEC_vs_Control_72h.csv \
        --out results/top10_common_degs.xlsx \
        --topn 10
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
    """Auto-detect column name from candidates."""
    m = {c.lower(): c for c in cols}
    for cand in candidates:
        if cand.lower() in m:
            return m[cand.lower()]
    if loose_keys:
        for c in cols:
            cl = c.lower()
            if all(k in cl for k in loose_keys):
                return c
    return None


def clean_subset(df, gene_col, logfc_col, p_col=None, padj_col=None):
    """Extract and clean relevant columns."""
    keep = [gene_col, logfc_col]
    if p_col:
        keep.append(p_col)
    if padj_col:
        keep.append(padj_col)
    d = df[keep].copy()
    d = d.dropna(subset=[gene_col, logfc_col])
    d[gene_col] = d[gene_col].astype(str).str.strip()
    d = d[d[gene_col] != ""]
    d[logfc_col] = pd.to_numeric(d[logfc_col], errors="coerce")
    d = d.dropna(subset=[logfc_col])
    return d


def dedup_keep_max_abs(df, gene_col, logfc_col):
    """Keep one row per gene with the largest |logFC|."""
    idx = df.groupby(gene_col)[logfc_col].apply(lambda s: s.abs().idxmax())
    return df.loc[idx]


def main():
    ap = argparse.ArgumentParser(
        description="Export top N common DEGs between 24h and 72h.")
    ap.add_argument("--deg24", required=True, type=Path)
    ap.add_argument("--deg72", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path, help="Output Excel file")
    ap.add_argument("--topn", type=int, default=10, help="Top N per direction (default: 10)")
    args = ap.parse_args()

    # Load data
    df24 = read_any(args.deg24)
    df72 = read_any(args.deg72)

    # Auto-detect columns
    gene_col = guess_col(df24.columns, ["Gene_Symbol", "SYMBOL", "Gene", "gene"], ["gene"])
    logfc_col = guess_col(df24.columns, ["logFC", "log2FC", "log2FoldChange"], ["log", "fc"])
    p_col = guess_col(df24.columns, ["P.Value", "pvalue", "p.value"], ["p", "val"])
    padj_col = guess_col(df24.columns, ["adj.P.Val", "padj", "FDR"], ["adj", "val"])

    if not gene_col or not logfc_col:
        raise SystemExit(
            f"Could not detect gene/logFC columns.\n"
            f"Detected: gene_col={gene_col}, logfc_col={logfc_col}")

    # Clean and rename
    d24 = clean_subset(df24, gene_col, logfc_col, p_col, padj_col).rename(columns={
        gene_col: "Gene", logfc_col: "logFC_24h",
        **(({p_col: "P.Value_24h"}) if p_col else {}),
        **(({padj_col: "adj.P.Val_24h"}) if padj_col else {}),
    })
    d72 = clean_subset(df72, gene_col, logfc_col, p_col, padj_col).rename(columns={
        gene_col: "Gene", logfc_col: "logFC_72h",
        **(({p_col: "P.Value_72h"}) if p_col else {}),
        **(({padj_col: "adj.P.Val_72h"}) if padj_col else {}),
    })

    # Deduplicate
    d24 = dedup_keep_max_abs(d24, "Gene", "logFC_24h")
    d72 = dedup_keep_max_abs(d72, "Gene", "logFC_72h")

    # Assign direction
    d24["dir_24h"] = np.where(d24["logFC_24h"] > 0, "up",
                        np.where(d24["logFC_24h"] < 0, "down", "zero"))
    d72["dir_72h"] = np.where(d72["logFC_72h"] > 0, "up",
                        np.where(d72["logFC_72h"] < 0, "down", "zero"))
    d24 = d24[d24["dir_24h"] != "zero"]
    d72 = d72[d72["dir_72h"] != "zero"]

    # Split by direction and find common
    d24u = d24[d24["dir_24h"] == "up"].set_index("Gene")
    d24d = d24[d24["dir_24h"] == "down"].set_index("Gene")
    d72u = d72[d72["dir_72h"] == "up"].set_index("Gene")
    d72d = d72[d72["dir_72h"] == "down"].set_index("Gene")

    # Join columns for common genes
    join_cols_24 = [c for c in ["logFC_24h", "P.Value_24h", "adj.P.Val_24h"] if c in d24u.columns]
    join_cols_72 = [c for c in ["logFC_72h", "P.Value_72h", "adj.P.Val_72h"] if c in d72u.columns]

    common_up = (d24u[join_cols_24]
                 .join(d72u[join_cols_72], how="inner")
                 .reset_index())
    common_down = (d24d[join_cols_24]
                   .join(d72d[join_cols_72], how="inner")
                   .reset_index())

    # Rank by mean logFC
    for df in (common_up, common_down):
        df["mean_logFC"] = df[["logFC_24h", "logFC_72h"]].mean(axis=1)

    top_up = common_up.sort_values("mean_logFC", ascending=False).head(args.topn)
    top_down = common_down.sort_values("mean_logFC", ascending=True).head(args.topn)

    # Reorder columns
    desired_order = ["Gene", "logFC_24h", "logFC_72h", "mean_logFC",
                     "P.Value_24h", "adj.P.Val_24h", "P.Value_72h", "adj.P.Val_72h"]
    top_up = top_up[[c for c in desired_order if c in top_up.columns]]
    top_down = top_down[[c for c in desired_order if c in top_down.columns]]

    # Write Excel
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(args.out, engine="xlsxwriter") as xw:
        top_up.to_excel(xw, index=False, sheet_name=f"Top{args.topn}_Common_Up")
        top_down.to_excel(xw, index=False, sheet_name=f"Top{args.topn}_Common_Down")

    print(f"Saved: {args.out}")
    print(f"  Common upregulated:   {len(common_up)} total, top {len(top_up)} exported")
    print(f"  Common downregulated: {len(common_down)} total, top {len(top_down)} exported")


if __name__ == "__main__":
    main()
