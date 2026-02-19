#!/usr/bin/env python3
"""Download GEO supplementary files for GSE249279.

Usage:
  python scripts/00_download_geo.py --config configs/config.yaml --skip-existing
"""
from __future__ import annotations

import argparse
from pathlib import Path

import requests
from tqdm import tqdm

from src.paper05.utils import read_yaml, ensure_dir

BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE249nnn/GSE249279/suppl"

def download_with_resume(url: str, dest: Path, chunk_size: int = 1024 * 1024) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    existing = dest.stat().st_size if dest.exists() else 0
    headers = {"Range": f"bytes={existing}-"} if existing > 0 else {}

    with requests.get(url, stream=True, headers=headers, timeout=120) as r:
        if r.status_code not in (200, 206):
            raise RuntimeError(f"Download failed: {url} status={r.status_code}")

        total = r.headers.get("Content-Length")
        total = int(total) + existing if total is not None else None

        mode = "ab" if existing > 0 else "wb"
        with open(dest, mode) as f, tqdm(
            total=total, initial=existing, unit="B", unit_scale=True, desc=dest.name
        ) as pbar:
            for chunk in r.iter_content(chunk_size=chunk_size):
                if chunk:
                    f.write(chunk)
                    pbar.update(len(chunk))

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="configs/config.yaml")
    ap.add_argument("--skip-existing", action="store_true")
    args = ap.parse_args()

    cfg = read_yaml(args.config)
    raw_dir = ensure_dir(cfg["paths"]["raw_dir"])
    files = cfg["download"]["files"]

    for fn in files:
        url = f"{BASE}/{fn}"
        dest = raw_dir / fn
        if args.skip_existing and dest.exists():
            print(f"[skip] {dest}")
            continue
        print(f"[get] {url}")
        download_with_resume(url, dest)

    print("[done] downloads complete")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
