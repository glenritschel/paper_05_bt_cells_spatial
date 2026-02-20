#!/usr/bin/env python3
"""
09_make_core_figures.py

Creates manuscript-ready core figures for Paper 05 (focused spatial replication):

Inputs
------
1) data/processed/spatial_programs/*.h5ad
   Must contain:
     - B_activation_score
     - Th2_score
     - Th17_score
     - total_counts
     - n_genes_by_counts
     - obsm['spatial'] (for optional spatial plots)
2) reports/b_th2_th17_correlations*.csv  (one or more files)

Outputs
-------
reports/figures/core/
  - scatter plots per section:
      {section}__Bact_vs_Th2__scatter.png
      {section}__Bact_vs_Th17__scatter.png
  - summary bar plots from correlation CSV:
      summary__Bact_vs_Th2__adj_spearman.png
      summary__Bact_vs_Th17__adj_spearman.png
  - optional spatial score panels for a chosen focus section:
      {focus}__spatial__B_activation_score.png
      {focus}__spatial__Th2_score.png
      {focus}__spatial__Th17_score.png

Notes
-----
- Uses Matplotlib only (no seaborn).
- Uses the same immune-only definition as 07_bt_th2_th17_correlations.py:
    immune spot = top-q in any of (B_activation, Th2, Th17)
- Performs "adjusted" scatter by residualizing both variables against
  total_counts and n_genes_by_counts.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
from scipy import stats

from src.paper05.utils import ensure_dir


def residualize(y: np.ndarray, X: np.ndarray) -> np.ndarray:
    """Return residuals of y ~ [1, X] least squares."""
    X_ = np.c_[np.ones(X.shape[0]), X]
    beta, *_ = np.linalg.lstsq(X_, y, rcond=None)
    return y - (X_ @ beta)


def immune_mask(b: np.ndarray, th2: np.ndarray, th17: np.ndarray, q: float) -> np.ndarray:
    b_thr = float(np.quantile(b, q))
    th2_thr = float(np.quantile(th2, q))
    th17_thr = float(np.quantile(th17, q))
    return (b >= b_thr) | (th2 >= th2_thr) | (th17 >= th17_thr)


def scatter_one(
    x: np.ndarray,
    y: np.ndarray,
    *,
    title: str,
    xlabel: str,
    ylabel: str,
    outpath: Path,
    annotate: str,
) -> None:
    plt.figure()
    plt.scatter(x, y, s=10, alpha=0.6)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    # Add a simple best-fit line (for visualization only)
    if len(x) >= 3:
        m, c = np.polyfit(x, y, 1)
        xs = np.linspace(np.min(x), np.max(x), 200)
        plt.plot(xs, m * xs + c)

    plt.text(
        0.02,
        0.98,
        annotate,
        transform=plt.gca().transAxes,
        ha="left",
        va="top",
        bbox=dict(boxstyle="round,pad=0.3", alpha=0.2),
    )

    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()


def plot_spatial_score(ad: sc.AnnData, score: str, outpath: Path) -> None:
    # use scanpy's spatial plot (squidpy version mismatch for img_key in some envs)
    sc.pl.spatial(ad, color=score, img_key="hires", size=1.3, show=False)
    fig = plt.gcf()
    fig.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)


def load_corr_csvs(paths: list[Path]) -> pd.DataFrame:
    dfs = []
    for p in paths:
        df = pd.read_csv(p)
        df["source_csv"] = p.name
        dfs.append(df)
    if not dfs:
        return pd.DataFrame()
    return pd.concat(dfs, axis=0, ignore_index=True)


def summary_barplot(df: pd.DataFrame, value_col: str, p_col: str, outpath: Path, title: str) -> None:
    # Expect sections A-D, but handle any
    d = df.sort_values("section").copy()
    sections = d["section"].astype(str).tolist()
    vals = d[value_col].astype(float).to_numpy()
    ps = d[p_col].astype(float).to_numpy()

    plt.figure()
    x = np.arange(len(sections))
    plt.bar(x, vals)
    plt.axhline(0.0, linewidth=1)

    plt.xticks(x, sections)
    plt.ylabel(value_col)
    plt.title(title)

    # p-value annotations
    for i, (v, p) in enumerate(zip(vals, ps)):
        plt.text(i, v, f"p={p:.3g}", ha="center", va="bottom" if v >= 0 else "top")

    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", default="data/processed/spatial_programs")
    ap.add_argument("--figdir", default="reports/figures/core")
    ap.add_argument("--corr-csv", default="reports/b_th2_th17_correlations.csv",
                    help="Path or glob, e.g. reports/b_th2_th17_correlations*.csv")
    ap.add_argument("--immune-only", action="store_true", help="Apply immune-only masking")
    ap.add_argument("--immune-q", type=float, default=0.9, help="Quantile for immune-only masking")
    ap.add_argument("--focus-section", default="D", help="Section letter to export spatial score maps for")
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()

    indir = Path(args.indir)
    figdir = ensure_dir(args.figdir)

    # correlation CSV(s)
    corr_paths = sorted(Path(".").glob(args.corr_csv)) if any(ch in args.corr_csv for ch in "*?[]") else [Path(args.corr_csv)]
    corr_paths = [p for p in corr_paths if p.exists()]
    corr_df = load_corr_csvs(corr_paths)

    # Load h5ads
    paths = sorted(indir.glob("*.h5ad"))
    if not paths:
        raise SystemExit(f"No inputs found: {indir}/*.h5ad")

    # Per-section scatterplots (adjusted residual scatter)
    for p in paths:
        section = p.stem
        ad = sc.read_h5ad(p)
        need = ["B_activation_score", "Th2_score", "Th17_score", "total_counts", "n_genes_by_counts"]
        missing = [k for k in need if k not in ad.obs.columns]
        if missing:
            raise SystemExit(f"{p.name} missing: {missing}")

        b = ad.obs["B_activation_score"].to_numpy()
        th2 = ad.obs["Th2_score"].to_numpy()
        th17 = ad.obs["Th17_score"].to_numpy()
        cov = ad.obs[["total_counts", "n_genes_by_counts"]].to_numpy()

        mask = np.ones(ad.n_obs, dtype=bool)
        if args.immune_only:
            mask = immune_mask(b, th2, th17, q=args.immune_q)

        b = b[mask]
        th2 = th2[mask]
        th17 = th17[mask]
        cov = cov[mask]

        # Adjusted residuals
        b_r = residualize(b, cov)
        th2_r = residualize(th2, cov)
        th17_r = residualize(th17, cov)

        # Stats for annotation (spearman on residuals)
        r_th2 = float(stats.spearmanr(b_r, th2_r).statistic)
        r_th17 = float(stats.spearmanr(b_r, th17_r).statistic)

        annotate_th2 = f"Spearman(adj) r={r_th2:.3f}\nN={len(b_r)}\nimmune_only={args.immune_only} q={args.immune_q}"
        annotate_th17 = f"Spearman(adj) r={r_th17:.3f}\nN={len(b_r)}\nimmune_only={args.immune_only} q={args.immune_q}"

        out1 = figdir / f"{section}__Bact_vs_Th2__scatter.png"
        out2 = figdir / f"{section}__Bact_vs_Th17__scatter.png"

        if args.overwrite or (not out1.exists()):
            scatter_one(
                b_r,
                th2_r,
                title=f"{section}: B activation vs Th2 (adjusted residuals)",
                xlabel="B_activation_score (residualized)",
                ylabel="Th2_score (residualized)",
                outpath=out1,
                annotate=annotate_th2,
            )
        if args.overwrite or (not out2.exists()):
            scatter_one(
                b_r,
                th17_r,
                title=f"{section}: B activation vs Th17 (adjusted residuals)",
                xlabel="B_activation_score (residualized)",
                ylabel="Th17_score (residualized)",
                outpath=out2,
                annotate=annotate_th17,
            )

        # Optional spatial panels for focus section
        if section == args.focus_section:
            for score in ["B_activation_score", "Th2_score", "Th17_score"]:
                sp_out = figdir / f"{section}__spatial__{score}.png"
                if args.overwrite or (not sp_out.exists()):
                    plot_spatial_score(ad, score, sp_out)

    # Summary barplots from correlation CSV (if present)
    if not corr_df.empty:
        # Prefer the "immune_only True" rows if present
        if "immune_only" in corr_df.columns:
            corr_df = corr_df.sort_values(["source_csv", "immune_only"], ascending=[True, False])

        # Choose one CSV to summarize: default the first one
        # If user passed multiple CSVs (q085/q090/q095), this can be extended later.
        src0 = corr_df["source_csv"].iloc[0]
        d0 = corr_df[corr_df["source_csv"] == src0].copy()

        # Barplots for adjusted spearman
        if "spearman_adj_B_vs_Th2" in d0.columns and "perm_p_adj_B_vs_Th2" in d0.columns:
            summary_barplot(
                d0,
                value_col="spearman_adj_B_vs_Th2",
                p_col="perm_p_adj_B_vs_Th2",
                outpath=figdir / "summary__Bact_vs_Th2__adj_spearman.png",
                title=f"Adjusted Spearman: B activation vs Th2 ({src0})",
            )

        if "spearman_adj_B_vs_Th17" in d0.columns and "perm_p_adj_B_vs_Th17" in d0.columns:
            summary_barplot(
                d0,
                value_col="spearman_adj_B_vs_Th17",
                p_col="perm_p_adj_B_vs_Th17",
                outpath=figdir / "summary__Bact_vs_Th17__adj_spearman.png",
                title=f"Adjusted Spearman: B activation vs Th17 ({src0})",
            )

    print(f"[done] wrote figures to {figdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

