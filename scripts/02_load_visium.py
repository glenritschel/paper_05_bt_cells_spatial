#!/usr/bin/env python3
"""Load Visium samples into AnnData using scanpy.read_visium()."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import scanpy as sc

from src.paper05.utils import read_yaml, ensure_dir

def find_sample_dirs(indir: Path) -> list[Path]:
    return sorted({p.parent for p in indir.rglob("spatial") if p.is_dir()})

def is_spaceranger_outs(d: Path) -> bool:
    return (d / "outs" / "spatial").is_dir() and (d / "outs" / "filtered_feature_bc_matrix.h5").exists()

def load_one(sample_dir: Path) -> sc.AnnData:
    vis_path = sample_dir / "outs" if is_spaceranger_outs(sample_dir) else sample_dir
    ad = sc.read_visium(path=str(vis_path))
    sc.pp.calculate_qc_metrics(ad, inplace=True)
    ad.uns["gse"] = "GSE249279"
    ad.uns["sample_dir"] = str(sample_dir)
    ad.uns["visium_path"] = str(vis_path)
    return ad

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="configs/config.yaml")
    ap.add_argument("--max-samples", type=int, default=0)
    args = ap.parse_args()

    cfg = read_yaml(args.config)
    interim_dir = Path(cfg["paths"]["interim_dir"])
    processed_dir = ensure_dir(Path(cfg["paths"]["processed_dir"]) / "spatial")

    sample_dirs = find_sample_dirs(interim_dir)
    if args.max_samples and args.max_samples > 0:
        sample_dirs = sample_dirs[: args.max_samples]

    out_recs = []
    for sd in sample_dirs:
        name = sd.name
        out = processed_dir / f"{name}.h5ad"
        print(f"[load] {sd} -> {out.name}")
        ad = load_one(sd)
        ad.write_h5ad(out)
        out_recs.append(
            {
                "sample_name": name,
                "sample_dir": str(sd),
                "n_spots": int(ad.n_obs),
                "n_genes": int(ad.n_vars),
                "out_h5ad": str(out),
            }
        )

    with open(processed_dir / "load_summary.json", "w") as f:
        json.dump(out_recs, f, indent=2)

    print(f"[done] Wrote {len(out_recs)} h5ad files to {processed_dir}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
