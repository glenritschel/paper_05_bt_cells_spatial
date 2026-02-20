#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

from src.paper05.utils import ensure_dir


def residualize(y: np.ndarray, X: np.ndarray) -> np.ndarray:
    X_ = np.c_[np.ones(X.shape[0]), X]
    beta, *_ = np.linalg.lstsq(X_, y, rcond=None)
    return y - (X_ @ beta)


def immune_mask(b: np.ndarray, th2: np.ndarray, th17: np.ndarray, q: float) -> np.ndarray:
    b_thr = float(np.quantile(b, q))
    th2_thr = float(np.quantile(th2, q))
    th17_thr = float(np.quantile(th17, q))
    return (b >= b_thr) | (th2 >= th2_thr) | (th17 >= th17_thr)


def to_dense_1d(x) -> np.ndarray:
    # x may be sparse matrix (n,1) or array
    if hasattr(x, "toarray"):
        x = x.toarray()
    x = np.asarray(x)
    return x.reshape(-1)


def nearest_distance(coords: np.ndarray, anchors: np.ndarray) -> np.ndarray:
    # O(n*m) is fine at Visium scale
    d = np.empty(coords.shape[0], dtype=float)
    for i in range(coords.shape[0]):
        diff = anchors - coords[i]
        d[i] = np.sqrt((diff * diff).sum(axis=1)).min()
    return d


def pick_anchors(lg: np.ndarray, q: float, topk: int, max_anchor_frac: float) -> tuple[np.ndarray, str]:
    """
    Choose anchor indices robustly.
    Start with quantile threshold.
    If too many anchors (degenerate threshold), fallback to top-k.
    """
    thr = float(np.quantile(lg, q))
    idx = np.where(lg >= thr)[0]

    if idx.size == 0:
        # fallback to top-k
        idx = np.argsort(lg)[-topk:]
        return idx, f"topk({topk})_fallback_empty_quantile"

    if idx.size / lg.size > max_anchor_frac:
        idx = np.argsort(lg)[-topk:]
        return idx, f"topk({topk})_fallback_too_many_quantile"

    return idx, f"quantile({q})"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", default="data/processed/spatial_programs")
    ap.add_argument("--outcsv", default="reports/lgals9_distance_effects_immune.csv")
    ap.add_argument("--lgals9-quantile", type=float, default=0.95)
    ap.add_argument("--topk", type=int, default=50)
    ap.add_argument("--max-anchor-frac", type=float, default=0.25, help="if quantile yields > this fraction anchors, fallback to topk")
    ap.add_argument("--immune-only", action="store_true")
    ap.add_argument("--immune-q", type=float, default=0.9)
    ap.add_argument("--adjust", action="store_true", help="residualize scores vs depth covariates before correlation")
    args = ap.parse_args()

    indir = Path(args.indir)
    paths = sorted(indir.glob("*.h5ad"))
    if not paths:
        raise SystemExit(f"No inputs found: {indir}/*.h5ad")

    ensure_dir(Path(args.outcsv).parent)

    rows = []
    for p in paths:
        section = p.stem
        ad = sc.read_h5ad(p)

        if "spatial" not in ad.obsm:
            raise SystemExit(f"{p.name} missing obsm['spatial']")
        if "T_suppression_exhaustion_score" not in ad.obs:
            raise SystemExit(f"{p.name} missing T_suppression_exhaustion_score")
        if "LGALS9" not in ad.var_names:
            print(f"[warn] {section}: LGALS9 not found; skipping")
            continue

        # immune-only filter (recommended for this question)
        mask = np.ones(ad.n_obs, dtype=bool)
        if args.immune_only:
            need = ["B_activation_score", "Th2_score", "Th17_score"]
            missing = [k for k in need if k not in ad.obs.columns]
            if missing:
                raise SystemExit(f"{p.name} missing for immune-only: {missing}")
            mask = immune_mask(
                ad.obs["B_activation_score"].to_numpy(),
                ad.obs["Th2_score"].to_numpy(),
                ad.obs["Th17_score"].to_numpy(),
                q=args.immune_q,
            )

        coords = ad.obsm["spatial"].astype(float)[mask]

        lg = to_dense_1d(ad[:, "LGALS9"].X)[mask]
        supp = ad.obs["T_suppression_exhaustion_score"].to_numpy()[mask]

        # Optional: adjust suppression vs depth (recommended)
        if args.adjust:
            cov = ad.obs[["total_counts", "n_genes_by_counts"]].to_numpy()[mask]
            supp = residualize(supp, cov)

        # anchor selection robust to degenerate quantiles
        anchor_idx_local, anchor_mode = pick_anchors(
            lg,
            q=args.lgals9_quantile,
            topk=args.topk,
            max_anchor_frac=args.max_anchor_frac,
        )
        if anchor_idx_local.size < 5:
            print(f"[warn] {section}: too few anchors ({anchor_idx_local.size}); skipping")
            continue

        anchors = coords[anchor_idx_local]
        dist = nearest_distance(coords, anchors)

        sp = float(stats.spearmanr(dist, supp).statistic)

        rows.append({
            "section": section,
            "n_spots_used": int(coords.shape[0]),
            "immune_only": bool(args.immune_only),
            "immune_q": float(args.immune_q),
            "lgals9_q": float(args.lgals9_quantile),
            "anchor_mode": anchor_mode,
            "n_anchors": int(anchor_idx_local.size),
            "spearman_dist_vs_Tsupp": sp,
            "adjusted": bool(args.adjust),
        })

    df = pd.DataFrame(rows).sort_values("section")
    df.to_csv(args.outcsv, index=False)
    print(f"[done] wrote {args.outcsv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

