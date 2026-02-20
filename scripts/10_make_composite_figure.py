#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc

from src.paper05.utils import ensure_dir


def to_dense_1d(x) -> np.ndarray:
    if hasattr(x, "toarray"):
        x = x.toarray()
    x = np.asarray(x)
    return x.reshape(-1)


def make_spatial_4panel(ad: sc.AnnData, outpath: Path) -> None:
    # Add LGALS9 expression to obs for plotting
    if "LGALS9" in ad.var_names:
        lg = to_dense_1d(ad[:, "LGALS9"].X)
        ad.obs["LGALS9_expr"] = lg
    else:
        ad.obs["LGALS9_expr"] = 0.0

    fig, axes = plt.subplots(2, 2, figsize=(10, 10))
    panels = [
        ("B_activation_score", "B activation"),
        ("Th2_score", "Th2 program"),
        ("Th17_score", "Th17 program"),
        ("LGALS9_expr", "LGALS9 expression"),
    ]

    for ax, (key, title) in zip(axes.ravel(), panels):
        sc.pl.spatial(ad, color=key, img_key="hires", size=1.3, show=False, ax=ax)
        ax.set_title(title)

    plt.tight_layout()
    fig.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)


def make_summary_barplot(df: pd.DataFrame, outpath: Path) -> None:
    # grouped bars: Th2 and Th17 adjusted spearman per section
    df = df.sort_values("section").copy()
    sections = df["section"].astype(str).tolist()

    th2 = df["spearman_adj_B_vs_Th2"].astype(float).to_numpy()
    th17 = df["spearman_adj_B_vs_Th17"].astype(float).to_numpy()

    p_th2 = df["perm_p_adj_B_vs_Th2"].astype(float).to_numpy()
    p_th17 = df["perm_p_adj_B_vs_Th17"].astype(float).to_numpy()

    x = np.arange(len(sections))
    w = 0.38

    fig = plt.figure(figsize=(10, 4))
    ax = plt.gca()

    ax.bar(x - w/2, th2, width=w, label="B vs Th2 (adj Spearman)")
    ax.bar(x + w/2, th17, width=w, label="B vs Th17 (adj Spearman)")
    ax.axhline(0.0, linewidth=1)

    ax.set_xticks(x)
    ax.set_xticklabels(sections)
    ax.set_ylabel("Adjusted Spearman r")
    ax.set_title("Immune-only (q=0.90): B activation vs Th2/Th17")
    ax.legend()

    # annotate p-values
    for i in range(len(sections)):
        ax.text(x[i]-w/2, th2[i], f"p={p_th2[i]:.3g}", ha="center",
                va="bottom" if th2[i] >= 0 else "top", fontsize=9)
        ax.text(x[i]+w/2, th17[i], f"p={p_th17[i]:.3g}", ha="center",
                va="bottom" if th17[i] >= 0 else "top", fontsize=9)

    plt.tight_layout()
    fig.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--spatial-programs", default="data/processed/spatial_programs")
    ap.add_argument("--corr-csv", default="reports/b_th2_th17_correlations.csv")
    ap.add_argument("--section", default="D")
    ap.add_argument("--outdir", default="reports/figures/composite")
    args = ap.parse_args()

    outdir = ensure_dir(args.outdir)

    # Load Section D
    h5 = Path(args.spatial_programs) / f"{args.section}.h5ad"
    if not h5.exists():
        raise SystemExit(f"Missing: {h5}")
    ad = sc.read_h5ad(h5)
    ad.var_names_make_unique()

    out_spatial = outdir / f"{args.section}__4panel__Bact_Th2_Th17_LGALS9.png"
    make_spatial_4panel(ad, out_spatial)

    # Summary bars
    df = pd.read_csv(args.corr_csv)
    out_bar = outdir / "summary__adj_spearman__Bact_vs_Th2_Th17.png"
    make_summary_barplot(df, out_bar)

    print(f"[done] wrote:\n  {out_spatial}\n  {out_bar}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

