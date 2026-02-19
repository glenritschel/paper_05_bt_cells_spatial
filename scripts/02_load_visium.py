#!/usr/bin/env python3
"""
Load Visium samples from extracted GSE249279_RAW into AnnData and write .h5ad outputs.

Heuristics:
- If sample_dir contains `outs/` (Space Ranger), point read_visium() at sample_dir/outs
- Else point read_visium() at sample_dir

Example:
  python scripts/02_load_visium.py --indir data/raw/extracted --outdir data/processed/spatial
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import scanpy as sc

def find_sample_dirs(indir: Path) -> list[Path]:
    # look for directories that have a spatial/ child
    dirs = set(p.parent for p in indir.rglob("spatial") if p.is_dir())
    return sorted(dirs)

def is_spaceranger_outs(d: Path) -> bool:
    return (d / "outs" / "spatial").is_dir() and (d / "outs" / "filtered_feature_bc_matrix.h5").exists()

def load_one(sample_dir: Path) -> sc.AnnData:
    if is_spaceranger_outs(sample_dir):
        vis_path = sample_dir / "outs"
    else:
        vis_path = sample_dir

    ad = sc.read_visium(
        path=str(vis_path),
        count_file="filtered_feature_bc_matrix.h5" if (vis_path / "filtered_feature_bc_matrix.h5").exists() else None
    )
    # Basic QC metrics in AnnData
    sc.pp.calculate_qc_metrics(ad, inplace=True)
    # keep provenance
    ad.uns["gse"] = "GSE249279"
    ad.uns["visium_path"] = str(vis_path)
    ad.uns["sample_dir"] = str(sample_dir)
    return ad

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", default="data/raw/extracted", help="Root dir containing extracted files")
    ap.add_argument("--outdir", default="data/processed/spatial", help="Output directory for .h5ad files")
    ap.add_argument("--max-samples", type=int, default=0, help="If >0, only load first N samples (debug)")
    args = ap.parse_args()

    indir = Path(args.indir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    sample_dirs = find_sample_dirs(indir)
    if args.max_samples and args.max_samples > 0:
        sample_dirs = sample_dirs[: args.max_samples]

    results = []
    print(f"[info] Found {len(sample_dirs)} candidate sample dirs (has spatial/)")
    for sd in sample_dirs:
        sample_name = sd.name
        out = outdir / f"{sample_name}.h5ad"
        print(f"[load] {sd} -> {out.name}")
        ad = load_one(sd)
        ad.write_h5ad(out)
        results.append({
            "sample_name": sample_name,
            "sample_dir": str(sd),
            "n_obs_spots": int(ad.n_obs),
            "n_vars_genes": int(ad.n_vars),
            "out_h5ad": str(out),
        })

    summary_path = outdir / "load_summary.json"
    with open(summary_path, "w") as f:
        json.dump(results, f, indent=2)

    print(f"[done] Wrote {len(results)} .h5ad files to {outdir}")
    print(f"[done] Summary: {summary_path}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
