#!/usr/bin/env python3
"""
Extract GSE249279_RAW.tar and verify Visium assets.

Writes:
- data/raw/manifest.json (by default)
- prints summary of detected Visium sample directories

Example:
  python scripts/01_extract_and_verify.py --raw-tar data/raw/GSE249279_RAW.tar --outdir data/raw/extracted
"""
from __future__ import annotations

import argparse
import json
import os
import tarfile
from pathlib import Path

VISIUM_MATRIX_HINTS = [
    "filtered_feature_bc_matrix.h5",
    "filtered_feature_bc_matrix",
    "raw_feature_bc_matrix.h5",
    "raw_feature_bc_matrix",
    "matrix.mtx.gz",
]

VISIUM_SPATIAL_HINTS = [
    "spatial/tissue_hires_image.png",
    "spatial/tissue_lowres_image.png",
    "spatial/scalefactors_json.json",
]

def safe_extract(tar: tarfile.TarFile, path: Path) -> None:
    # basic path traversal guard
    for member in tar.getmembers():
        member_path = path / member.name
        if not str(member_path.resolve()).startswith(str(path.resolve())):
            raise RuntimeError(f"Unsafe path in tar: {member.name}")
    tar.extractall(path)

def build_inventory(root: Path) -> dict:
    files = []
    for p in root.rglob("*"):
        if p.is_file():
            files.append(str(p.relative_to(root)))
    files.sort()
    return {"root": str(root), "n_files": len(files), "files": files}

def detect_visium_samples(root: Path) -> list[dict]:
    """
    Heuristic: a Visium sample dir contains:
      - a matrix asset (h5 or mtx)
      - a spatial/ folder with expected files
    """
    candidates = []
    for spatial_dir in root.rglob("spatial"):
        if not spatial_dir.is_dir():
            continue
        sample_dir = spatial_dir.parent
        # look for matrix hints in sample_dir (or typical subpaths)
        has_matrix = False
        for hint in VISIUM_MATRIX_HINTS:
            if (sample_dir / hint).exists():
                has_matrix = True
                break
        # also check common Space Ranger layout: sample_dir/outs/filtered_feature_bc_matrix.h5
        if not has_matrix and (sample_dir / "outs" / "filtered_feature_bc_matrix.h5").exists():
            has_matrix = True

        has_spatial = any((sample_dir / h).exists() for h in VISIUM_SPATIAL_HINTS)
        # tissue_positions list can be in multiple formats
        has_positions = any((sample_dir / "spatial" / fn).exists() for fn in [
            "tissue_positions_list.csv",
            "tissue_positions.csv",
            "tissue_positions.parquet",
        ])

        if has_matrix and (has_spatial or has_positions):
            candidates.append({
                "sample_dir": str(sample_dir),
                "has_matrix": has_matrix,
                "has_spatial_images_and_scalefactors": has_spatial,
                "has_tissue_positions": has_positions,
                "matrix_paths_found": [
                    str((sample_dir / hint).relative_to(root))
                    for hint in VISIUM_MATRIX_HINTS
                    if (sample_dir / hint).exists()
                ] + ([
                    str((sample_dir / "outs" / "filtered_feature_bc_matrix.h5").relative_to(root))
                ] if (sample_dir / "outs" / "filtered_feature_bc_matrix.h5").exists() else []),
            })
    # de-dup by sample_dir
    uniq = {}
    for c in candidates:
        uniq[c["sample_dir"]] = c
    return sorted(uniq.values(), key=lambda x: x["sample_dir"])

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-tar", required=True, help="Path to GSE249279_RAW.tar")
    ap.add_argument("--outdir", default="data/raw/extracted", help="Directory to extract into")
    ap.add_argument("--manifest", default="data/raw/manifest.json", help="Where to write manifest JSON")
    ap.add_argument("--skip-extract", action="store_true", help="Assume already extracted; just verify")
    args = ap.parse_args()

    raw_tar = Path(args.raw_tar)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if not args.skip_extract:
        print(f"[info] Extracting: {raw_tar} -> {outdir}")
        with tarfile.open(raw_tar, "r:*") as tar:
            safe_extract(tar, outdir)
        print("[done] extract complete")
    else:
        print("[info] Skipping extract")

    inv = build_inventory(outdir)
    vis = detect_visium_samples(outdir)

    manifest = {
        "dataset": "GSE249279",
        "raw_tar": str(raw_tar),
        "extracted_dir": str(outdir),
        "inventory": {"n_files": inv["n_files"]},
        "visium_samples_detected": vis,
    }

    Path(args.manifest).parent.mkdir(parents=True, exist_ok=True)
    with open(args.manifest, "w") as f:
        json.dump(manifest, f, indent=2)

    print(f"[info] Inventory: {inv['n_files']} files under {outdir}")
    print(f"[info] Detected Visium samples: {len(vis)}")
    for v in vis:
        print(f"  - {v['sample_dir']}")
        if v["matrix_paths_found"]:
            print(f"      matrix: {v['matrix_paths_found'][0]}")
        print(f"      spatial_images_scalefactors={v['has_spatial_images_and_scalefactors']} positions={v['has_tissue_positions']}")
    print(f"[done] Wrote manifest: {args.manifest}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
