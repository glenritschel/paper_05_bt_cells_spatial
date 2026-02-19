#!/usr/bin/env python3
"""Extract GSE249279_RAW.tar and detect Visium sample directories."""
from __future__ import annotations

import argparse
import tarfile
from pathlib import Path
from typing import Any, Dict, List

from src.paper05.utils import read_yaml, ensure_dir, write_json

VISIUM_MATRIX_HINTS = ["filtered_feature_bc_matrix.h5", "raw_feature_bc_matrix.h5", "matrix.mtx.gz"]

def safe_extract(tar: tarfile.TarFile, path: Path) -> None:
    for member in tar.getmembers():
        member_path = path / member.name
        if not str(member_path.resolve()).startswith(str(path.resolve())):
            raise RuntimeError(f"Unsafe path in tar: {member.name}")
    tar.extractall(path)

def detect_visium_samples(root: Path) -> List[Dict[str, Any]]:
    samples: List[Dict[str, Any]] = []
    for spatial_dir in root.rglob("spatial"):
        if not spatial_dir.is_dir():
            continue
        sample_dir = spatial_dir.parent

        matrix_paths: List[str] = []
        for hint in VISIUM_MATRIX_HINTS:
            if (sample_dir / hint).exists():
                matrix_paths.append(str(sample_dir / hint))
        if (sample_dir / "outs" / "filtered_feature_bc_matrix.h5").exists():
            matrix_paths.append(str(sample_dir / "outs" / "filtered_feature_bc_matrix.h5"))

        has_scalefactors = (sample_dir / "spatial" / "scalefactors_json.json").exists()
        has_hires = (sample_dir / "spatial" / "tissue_hires_image.png").exists()
        has_positions = any(
            (sample_dir / "spatial" / fn).exists()
            for fn in ["tissue_positions_list.csv", "tissue_positions.csv", "tissue_positions.parquet"]
        )

        if matrix_paths and (has_scalefactors or has_positions or has_hires):
            samples.append(
                {
                    "sample_dir": str(sample_dir),
                    "matrix_paths": matrix_paths,
                    "has_scalefactors": has_scalefactors,
                    "has_hires_image": has_hires,
                    "has_tissue_positions": has_positions,
                }
            )

    uniq = {s["sample_dir"]: s for s in samples}
    return sorted(uniq.values(), key=lambda x: x["sample_dir"])

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="configs/config.yaml")
    ap.add_argument("--skip-extract", action="store_true")
    args = ap.parse_args()

    cfg = read_yaml(args.config)
    raw_dir = Path(cfg["paths"]["raw_dir"])
    interim_dir = ensure_dir(cfg["paths"]["interim_dir"])
    tar_path = raw_dir / "GSE249279_RAW.tar"

    if not args.skip_extract:
        print(f"[info] Extracting {tar_path} -> {interim_dir}")
        with tarfile.open(tar_path, "r:*") as tar:
            safe_extract(tar, interim_dir)

    n_files = sum(1 for p in interim_dir.rglob("*") if p.is_file())
    vis = detect_visium_samples(interim_dir)

    manifest = {
        "dataset": cfg["dataset"]["gse"],
        "raw_tar": str(tar_path),
        "interim_dir": str(interim_dir),
        "n_files": n_files,
        "visium_samples": vis,
    }
    write_json(manifest, raw_dir / "manifest.json")

    print(f"[info] Extracted file count: {n_files}")
    print(f"[info] Detected Visium samples: {len(vis)}")
    for v in vis:
        print(f"  - {v['sample_dir']}")
    print(f"[done] Wrote {raw_dir/'manifest.json'}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
