#!/usr/bin/env python3
"""
Score B/T/Plasma marker gene sets on preprocessed Visium sections and export spatial plots.

Inputs:
  data/processed/spatial_prepped/*.h5ad   (preferred)
  (fallback) data/processed/spatial/*.h5ad

Outputs:
  - writes scores into ad.obs
  - saves updated .h5ad alongside input (or to --outdir)
  - exports PNGs to reports/figures/bt_scores/

Usage:
  python scripts/04_score_bt_markers.py --config configs/config.yaml
"""
from __future__ import annotations

import argparse
from pathlib import Path

import scanpy as sc
import matplotlib.pyplot as plt

from src.paper05.utils import read_yaml, ensure_dir


# Conservative, widely-used markers (you can expand later)
DEFAULT_MARKERS = {
    "B_score": ["MS4A1", "CD79A", "CD79B", "CD74", "HLA-DRA"],
    "T_score": ["CD3D", "CD3E", "CD3G", "TRAC", "LCK"],
    "Plasma_score": ["MZB1", "XBP1", "SDC1", "JCHAIN", "IGHG1"],
    # optional "sanity" controls
    "Fibro_score": ["COL1A1", "COL1A2", "DCN", "LUM"],
    "Endo_score": ["KDR", "PECAM1", "VWF", "EMCN"],
    "Kerat_score": ["KRT14", "KRT5", "KRT1", "KRT10"],
}


def score_sets(ad: sc.AnnData, markers: dict[str, list[str]]) -> None:
    # score_genes ignores missing genes; we also log how many are present
    present = set(ad.var_names)
    for score_name, genes in markers.items():
        genes_present = [g for g in genes if g in present]
        ad.uns.setdefault("marker_scoring", {})[score_name] = {
            "requested": genes,
            "present": genes_present,
            "n_present": len(genes_present),
        }
        sc.tl.score_genes(ad, gene_list=genes_present, score_name=score_name, use_raw=False)


def plot_spatial_scores(ad: sc.AnnData, outdir: Path, section_name: str) -> None:
    outdir.mkdir(parents=True, exist_ok=True)

    # Core B/T/Plasma
    core = ["B_score", "T_score", "Plasma_score"]
    for k in core:
        if k in ad.obs:
            sc.pl.spatial(
                ad,
                color=k,
                img_key="hires",
                size=1.3,
                show=False,
            )
            plt.tight_layout()
            plt.savefig(outdir / f"{section_name}__{k}.png", dpi=200)
            plt.close()

    # Side-by-side (B vs T)
    if "B_score" in ad.obs and "T_score" in ad.obs:
        sc.pl.spatial(
            ad,
            color=["B_score", "T_score"],
            img_key="hires",
            size=1.3,
            show=False,
        )
        plt.tight_layout()
        plt.savefig(outdir / f"{section_name}__B_vs_T.png", dpi=200)
        plt.close()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="configs/config.yaml")
    ap.add_argument("--indir", default="data/processed/spatial_prepped")
    ap.add_argument("--fallback-indir", default="data/processed/spatial")
    ap.add_argument("--outdir", default="data/processed/spatial_scored")
    ap.add_argument("--figdir", default="reports/figures/bt_scores")
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()

    _cfg = read_yaml(args.config)  # currently unused, kept for consistency

    indir = Path(args.indir)
    inputs = sorted(indir.glob("*.h5ad"))
    if not inputs:
        indir = Path(args.fallback_indir)
        inputs = sorted(indir.glob("*.h5ad"))

    if not inputs:
        raise SystemExit(f"No inputs found in {args.indir} or {args.fallback_indir}")

    outdir = ensure_dir(args.outdir)
    figdir = ensure_dir(args.figdir)

    for p in inputs:
        section = p.stem
        out = outdir / p.name
        if out.exists() and not args.overwrite:
            print(f"[skip] {out}")
            continue

        print(f"[score] {p} -> {out}")
        ad = sc.read_h5ad(p)
        ad.var_names_make_unique()

        # Ensure log1p exists for score_genes (score_genes assumes expression scale; log1p is OK)
        # If you prefer raw counts scoring, we can store counts in layers and score on that later.
        score_sets(ad, DEFAULT_MARKERS)

        # Save scored object
        ad.write_h5ad(out)

        # Export spatial plots
        plot_spatial_scores(ad, figdir, section)

        # Print presence summary (useful for debugging gene naming conventions)
        ms = ad.uns.get("marker_scoring", {})
        core = {k: ms.get(k, {}).get("n_present", 0) for k in ["B_score", "T_score", "Plasma_score"]}
        print(f"  [markers present] {core}")

    print(f"[done] Wrote scored h5ad files to {outdir} and figures to {figdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

