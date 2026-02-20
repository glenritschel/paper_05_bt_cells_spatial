#!/usr/bin/env python3
"""
Preprocess Visium sections (.h5ad) for analysis-ready spatial work.

Inputs:
  data/processed/spatial/*.h5ad

Outputs:
  data/processed/spatial_prepped/*.h5ad

Steps (per section):
  - basic spot filters (min_counts, min_genes)
  - normalize_total + log1p
  - HVGs (optional; stored in ad.var['highly_variable'])
  - PCA / neighbors / UMAP
  - Leiden clustering

Usage:
  python scripts/03_preprocess_visium.py --config configs/config.yaml
"""
from __future__ import annotations

import argparse
from pathlib import Path

import scanpy as sc

from src.paper05.utils import read_yaml, ensure_dir


def preprocess_one(
    ad: sc.AnnData,
    *,
    min_counts: int,
    min_genes: int,
    target_sum: float,
    n_top_hvgs: int,
    n_pcs: int,
    n_neighbors: int,
    leiden_resolution: float,
) -> sc.AnnData:
    # Basic filtering
    sc.pp.filter_cells(ad, min_counts=min_counts)
    sc.pp.filter_cells(ad, min_genes=min_genes)

    # Normalize + log
    sc.pp.normalize_total(ad, target_sum=target_sum)
    sc.pp.log1p(ad)

    # HVGs (useful for PCA/UMAP/clustering). Keep all genes but mark HVGs.
    if n_top_hvgs and n_top_hvgs > 0:
        # seurat_v3 can require raw counts; after log1p it still usually works, but we keep it simple:
        sc.pp.highly_variable_genes(ad, flavor="seurat", n_top_genes=n_top_hvgs, subset=False)

    # PCA / neighbors / UMAP
    sc.tl.pca(ad, n_comps=n_pcs, use_highly_variable=bool("highly_variable" in ad.var))
    sc.pp.neighbors(ad, n_neighbors=n_neighbors, n_pcs=n_pcs)
    sc.tl.umap(ad)

    # Clustering
    sc.tl.leiden(ad, resolution=leiden_resolution, key_added="leiden")

    return ad


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="configs/config.yaml")
    ap.add_argument("--indir", default="data/processed/spatial")
    ap.add_argument("--outdir", default="data/processed/spatial_prepped")
    ap.add_argument("--min-counts", type=int, default=500)
    ap.add_argument("--min-genes", type=int, default=200)
    ap.add_argument("--target-sum", type=float, default=1e4)
    ap.add_argument("--n-top-hvgs", type=int, default=3000)
    ap.add_argument("--n-pcs", type=int, default=30)
    ap.add_argument("--n-neighbors", type=int, default=15)
    ap.add_argument("--leiden-resolution", type=float, default=0.5)
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()

    _cfg = read_yaml(args.config)  # currently unused, but kept for consistency

    indir = Path(args.indir)
    outdir = ensure_dir(args.outdir)

    inputs = sorted(indir.glob("*.h5ad"))
    if not inputs:
        raise SystemExit(f"No inputs found: {indir}/*.h5ad")

    for p in inputs:
        out = outdir / p.name
        if out.exists() and not args.overwrite:
            print(f"[skip] {out}")
            continue

        print(f"[prep] {p} -> {out}")
        ad = sc.read_h5ad(p)

        # Ensure duplicates are fixed (should already be true, but safe)
        ad.var_names_make_unique()

        # Compute QC metrics if missing
        if "total_counts" not in ad.obs or "n_genes_by_counts" not in ad.obs:
            sc.pp.calculate_qc_metrics(ad, inplace=True)

        ad = preprocess_one(
            ad,
            min_counts=args.min_counts,
            min_genes=args.min_genes,
            target_sum=args.target_sum,
            n_top_hvgs=args.n_top_hvgs,
            n_pcs=args.n_pcs,
            n_neighbors=args.n_neighbors,
            leiden_resolution=args.leiden_resolution,
        )

        ad.uns["preprocess"] = {
            "min_counts": args.min_counts,
            "min_genes": args.min_genes,
            "target_sum": args.target_sum,
            "n_top_hvgs": args.n_top_hvgs,
            "n_pcs": args.n_pcs,
            "n_neighbors": args.n_neighbors,
            "leiden_resolution": args.leiden_resolution,
        }

        ad.write_h5ad(out)

    print(f"[done] Wrote prepped h5ad files to {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

