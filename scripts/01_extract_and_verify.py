#!/usr/bin/env python3
"""Extract GSE249279_RAW.tar and detect Visium sample directories."""
from __future__ import annotations

import argparse
import tarfile
import hashlib
import time
import json

from pathlib import Path
from typing import Any, Dict, List

from src.paper05.utils import read_yaml, ensure_dir, write_json

VISIUM_MATRIX_HINTS = ["filtered_feature_bc_matrix.h5", "raw_feature_bc_matrix.h5", "matrix.mtx.gz"]


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        while True:
            b = f.read(chunk_size)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def stamp_write(stamp: Path, payload: dict) -> None:
    stamp.parent.mkdir(parents=True, exist_ok=True)
    stamp.write_text(json.dumps(payload, indent=2) + "\n")


def stamp_ok(stamp: Path, expected: dict) -> bool:
    if not stamp.exists():
        return False
    try:
        got = json.loads(stamp.read_text())
    except Exception:
        return False
    return got == expected


def safe_extract(tar: tarfile.TarFile, path: Path) -> None:
    for member in tar.getmembers():
        member_path = path / member.name
        if not str(member_path.resolve()).startswith(str(path.resolve())):
            raise RuntimeError(f"Unsafe path in tar: {member.name}")
    tar.extractall(path)


def extract_nested_archives(root: Path, force: bool = False) -> int:
    """
    Extract nested .tar, .tar.gz, .tgz files inside `root`
    into sibling directories named <archive>.extracted

    Returns number of archives extracted.
    """
    import tarfile

    n_extracted = 0

    for archive in sorted(root.rglob("*")):
        if not archive.is_file():
            continue

        name = archive.name.lower()
        if not (name.endswith(".tar") or name.endswith(".tar.gz") or name.endswith(".tgz")):
            continue

        outdir = archive.parent / f"{archive.name}.extracted"

        # Skip if already extracted and non-empty
        if outdir.exists() and any(outdir.rglob("*")):
            continue

        print(f"[info] Extracting nested archive {archive} -> {outdir}")
        outdir.mkdir(parents=True, exist_ok=True)

        arch_sha = sha256_file(archive)
        stamp = outdir / ".extract.stamp.json"
        expected = {"archive": str(archive), "sha256": arch_sha}

        if (not force) and stamp_ok(stamp, expected):
            continue

        # extract...
        stamp_write(stamp, expected)

        with tarfile.open(archive, "r:*") as tar:
            safe_extract(tar, outdir)

        n_extracted += 1

    return n_extracted


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
    ap.add_argument("--force", action="store_true", help="Force re-extraction ignoring stamps")
    args = ap.parse_args()

    cfg = read_yaml(args.config)
    raw_dir = Path(cfg["paths"]["raw_dir"])
    interim_dir = ensure_dir(cfg["paths"]["interim_dir"])
    tar_path = raw_dir / "GSE249279_RAW.tar"

    raw_sha = sha256_file(tar_path)
    outer_stamp = interim_dir / ".outer_extract.stamp.json"
    outer_expected = {"tar": str(tar_path), "sha256": raw_sha}

    if not args.skip_extract:
        if (not args.force) and stamp_ok(outer_stamp, outer_expected):
            print("[info] Outer tar already extracted (stamp match); skipping.")
        else:
            print(f"[info] Extracting {tar_path} -> {interim_dir}")
            # optional: clear interim_dir except .gitkeep and extracted dirs
            with tarfile.open(tar_path, "r:*") as tar:
                safe_extract(tar, interim_dir)
            stamp_write(outer_stamp, outer_expected)

    nested = extract_nested_archives(interim_dir, force=args.force)
    print(f"[info] Nested archives extracted: {nested}")

    print(f"[info] Nested archives extracted: {nested}")

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
