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
    # y ~ X via least squares; return residuals
    X_ = np.c_[np.ones(X.shape[0]), X]
    beta, *_ = np.linalg.lstsq(X_, y, rcond=None)
    yhat = X_ @ beta
    return y - yhat


def perm_pvalue(x, y, n=5000, seed=0):
    rng = np.random.default_rng(seed)
    obs = stats.spearmanr(x, y).statistic
    cnt = 0
    for _ in range(n):
        yp = rng.permutation(y)
        r = stats.spearmanr(x, yp).statistic
        if abs(r) >= abs(obs):
            cnt += 1
    return obs, (cnt + 1) / (n + 1)

def section_metrics(
    ad: sc.AnnData,
    section: str,
    perms: int,
    immune_only: bool = False,
    immune_q: float = 0.9,
) -> dict:
    b = ad.obs["B_activation_score"].to_numpy()
    th2 = ad.obs["Th2_score"].to_numpy()
    th17 = ad.obs["Th17_score"].to_numpy()

    # covariates
    cov = ad.obs[["total_counts", "n_genes_by_counts"]].to_numpy()

    mask = np.ones(ad.n_obs, dtype=bool)
    if immune_only:
        b_thr = float(np.quantile(b, immune_q))
        th2_thr = float(np.quantile(th2, immune_q))
        th17_thr = float(np.quantile(th17, immune_q))
        mask = (b >= b_thr) | (th2 >= th2_thr) | (th17 >= th17_thr)

    b = b[mask]
    th2 = th2[mask]
    th17 = th17[mask]
    cov = cov[mask]

    # residualize to reduce depth effects
    b_r = residualize(b, cov)
    th2_r = residualize(th2, cov)
    th17_r = residualize(th17, cov)

    out = {"section": section, "n_spots": int(b.shape[0]), "immune_only": bool(immune_only), "immune_q": float(immune_q)}

    # raw correlations
    out["spearman_B_vs_Th2"] = float(stats.spearmanr(b, th2).statistic)
    out["spearman_B_vs_Th17"] = float(stats.spearmanr(b, th17).statistic)
    out["pearson_B_vs_Th2"] = float(np.corrcoef(b, th2)[0, 1])
    out["pearson_B_vs_Th17"] = float(np.corrcoef(b, th17)[0, 1])

    # adjusted correlations (spearman on residuals)
    out["spearman_adj_B_vs_Th2"] = float(stats.spearmanr(b_r, th2_r).statistic)
    out["spearman_adj_B_vs_Th17"] = float(stats.spearmanr(b_r, th17_r).statistic)

    # permutation p-values on adjusted (more conservative)
    r, p = perm_pvalue(b_r, th2_r, n=perms, seed=1)
    out["perm_spearman_adj_B_vs_Th2"] = float(r)
    out["perm_p_adj_B_vs_Th2"] = float(p)

    r, p = perm_pvalue(b_r, th17_r, n=perms, seed=2)
    out["perm_spearman_adj_B_vs_Th17"] = float(r)
    out["perm_p_adj_B_vs_Th17"] = float(p)

    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", default="data/processed/spatial_programs")
    ap.add_argument("--outcsv", default="reports/b_th2_th17_correlations.csv")
    ap.add_argument("--perms", type=int, default=2000)
    ap.add_argument("--immune-only", action="store_true")
    ap.add_argument("--immune-q", type=float, default=0.9)

    args = ap.parse_args()

    indir = Path(args.indir)
    paths = sorted(indir.glob("*.h5ad"))
    if not paths:
        raise SystemExit(f"No inputs found: {indir}/*.h5ad")

    ensure_dir(Path(args.outcsv).parent)

    rows = []
    for p in paths:
        ad = sc.read_h5ad(p)
        need = ["B_activation_score", "Th2_score", "Th17_score", "total_counts", "n_genes_by_counts"]
        missing = [k for k in need if k not in ad.obs.columns]
        if missing:
            raise SystemExit(f"{p.name} missing: {missing}")

        rows.append(section_metrics(
            ad,
            p.stem,
            perms=args.perms,
            immune_only=args.immune_only,
            immune_q=args.immune_q
        ))

    df = pd.DataFrame(rows).sort_values("section")
    df.to_csv(args.outcsv, index=False)
    print(f"[done] wrote {args.outcsv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

