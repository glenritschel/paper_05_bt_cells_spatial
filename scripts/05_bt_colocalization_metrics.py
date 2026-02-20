#!/usr/bin/env python3
"""
Compute B<->T colocalization metrics per Visium section.

Inputs:
  data/processed/spatial_scored/*.h5ad  (expects B_score and T_score in ad.obs)

Outputs:
  reports/bt_metrics.csv
  reports/figures/bt_metrics/*

Metrics per section:
  - Pearson/Spearman corr(B_score, T_score)
  - overlap of top decile (B_high & T_high)
  - Squidpy neighborhood enrichment on bt_bin (B_high / T_high / Other)
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt

from src.paper05.utils import read_yaml, ensure_dir


def assign_bins(ad: sc.AnnData, q: float = 0.9) -> None:
    bq = float(np.quantile(ad.obs["B_score"].values, q))
    tq = float(np.quantile(ad.obs["T_score"].values, q))
    b_high = ad.obs["B_score"].values >= bq
    t_high = ad.obs["T_score"].values >= tq

    bt = np.full(ad.n_obs, "Other", dtype=object)
    bt[b_high & ~t_high] = "B_high"
    bt[~b_high & t_high] = "T_high"
    bt[b_high & t_high] = "B_and_T_high"

    ad.obs["bt_bin"] = pd.Categorical(bt, categories=["Other", "B_high", "T_high", "B_and_T_high"])


def corr_metrics(ad: sc.AnnData) -> dict:
    b = pd.Series(ad.obs["B_score"].values)
    t = pd.Series(ad.obs["T_score"].values)
    return {
        "pearson_r": float(b.corr(t, method="pearson")),
        "spearman_r": float(b.corr(t, method="spearman")),
    }


def top_decile_overlap(ad: sc.AnnData, q: float = 0.9) -> dict:
    bq = float(np.quantile(ad.obs["B_score"].values, q))
    tq = float(np.quantile(ad.obs["T_score"].values, q))
    b_high = ad.obs["B_score"].values >= bq
    t_high = ad.obs["T_score"].values >= tq

    n = ad.n_obs
    nb = int(b_high.sum())
    nt = int(t_high.sum())
    nbt = int((b_high & t_high).sum())

    # conditional fractions
    frac_b_also_t = float(nbt / nb) if nb else np.nan
    frac_t_also_b = float(nbt / nt) if nt else np.nan
    frac_both = float(nbt / n) if n else np.nan

    return {
        "q": q,
        "n_spots": int(n),
        "n_B_high": nb,
        "n_T_high": nt,
        "n_B_and_T_high": nbt,
        "frac_B_high_also_T_high": frac_b_also_t,
        "frac_T_high_also_B_high": frac_t_also_b,
        "frac_spots_B_and_T_high": frac_both,
    }


def plot_bt_bins(ad: sc.AnnData, figdir: Path, section: str) -> None:
    sc.pl.spatial(ad, color="bt_bin", img_key="hires", size=1.3, show=False)
    fig = plt.gcf()
    fig.savefig(figdir / f"{section}__bt_bin.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def run_nhood_enrichment(ad: sc.AnnData, figdir: Path, section: str) -> dict:
    # Build spatial neighbors and neighborhood enrichment.
    # coord_type="grid" is appropriate for Visium spot lattice.
    sq.gr.spatial_neighbors(ad, coord_type="grid")
    sq.gr.nhood_enrichment(ad, cluster_key="bt_bin")

    # Plot enrichment heatmap
    sq.pl.nhood_enrichment(ad, cluster_key="bt_bin", figsize=(4, 4), show=False)

    # Some squidpy versions return None; use current figure instead.
    fig = plt.gcf()
    fig.savefig(figdir / f"{section}__nhood_enrichment.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    # Extract z-score matrix (stored in ad.uns)
    z = ad.uns.get("bt_bin_nhood_enrichment", {}).get("zscore", None)
    if z is None:
        return {}

    # z is a square matrix aligned to categories
    cats = list(ad.obs["bt_bin"].cat.categories)
    zdf = pd.DataFrame(z, index=cats, columns=cats)

    # Provide a couple of “headline” pairs
    out = {}
    for a, b in [("B_high", "T_high"), ("B_high", "B_high"), ("T_high", "T_high"), ("B_and_T_high", "Other")]:
        if a in zdf.index and b in zdf.columns:
            out[f"z_{a}__adj__{b}"] = float(zdf.loc[a, b])
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="configs/config.yaml")
    ap.add_argument("--indir", default="data/processed/spatial_scored")
    ap.add_argument("--outcsv", default="reports/bt_metrics.csv")
    ap.add_argument("--figdir", default="reports/figures/bt_metrics")
    ap.add_argument("--q", type=float, default=0.9, help="quantile for high-score binning")
    ap.add_argument("--overwrite-figs", action="store_true")
    args = ap.parse_args()

    _cfg = read_yaml(args.config)  # reserved for future use

    indir = Path(args.indir)
    paths = sorted(indir.glob("*.h5ad"))
    if not paths:
        raise SystemExit(f"No inputs found: {indir}/*.h5ad")

    figdir = ensure_dir(args.figdir)
    outcsv = Path(args.outcsv)
    outcsv.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    for p in paths:
        section = p.stem
        print(f"[metrics] {p.name}")
        ad = sc.read_h5ad(p)

        if "B_score" not in ad.obs or "T_score" not in ad.obs:
            raise SystemExit(f"{p} missing B_score/T_score. Run score step first.")

        assign_bins(ad, q=args.q)

        row = {"section": section}
        row.update(corr_metrics(ad))
        row.update(top_decile_overlap(ad, q=args.q))

        # figures
        bt_fig = figdir / f"{section}__bt_bin.png"
        nh_fig = figdir / f"{section}__nhood_enrichment.png"
        if args.overwrite_figs or (not bt_fig.exists()):
            plot_bt_bins(ad, figdir, section)

        # neighborhood enrichment (computational; run each time)
        row.update(run_nhood_enrichment(ad, figdir, section))

        rows.append(row)

    df = pd.DataFrame(rows).sort_values("section")
    df.to_csv(outcsv, index=False)
    print(f"[done] wrote {outcsv} and figures in {figdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

