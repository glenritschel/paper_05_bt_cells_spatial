#!/usr/bin/env python3
"""
Download GEO series supplementary files for GSE249279.

Supports:
- direct NCBI FTP/HTTPS download (default)
- GEOparse-driven discovery (optional)

Examples:
  python scripts/00_download_geo.py --outdir data/raw
  python scripts/00_download_geo.py --method geoparse --outdir data/raw
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path
import requests
from tqdm import tqdm

GSE = "GSE249279"

# NCBI GEO series structure: .../geo/series/GSE249nnn/GSE249279/suppl/<file>
BASE_FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE249nnn/GSE249279/suppl"

FILES = [
    "GSE249279_RAW.tar",
    "GSE249279_scRNA-seq_merged_raw_counts.h5",
]

def download_with_resume(url: str, dest: Path, chunk_size: int = 1024 * 1024) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    existing = dest.stat().st_size if dest.exists() else 0

    headers = {}
    if existing > 0:
        headers["Range"] = f"bytes={existing}-"

    with requests.get(url, stream=True, headers=headers, timeout=60) as r:
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

def direct_urls() -> dict[str, str]:
    return {fn: f"{BASE_FTP}/{fn}" for fn in FILES}

def geoparse_urls() -> dict[str, str]:
    try:
        import GEOparse  # type: ignore
    except Exception as e:
        raise RuntimeError("GEOparse not installed. Install from requirements.txt or use --method direct.") from e

    gse = GEOparse.get_GEO(GSE, destdir=".")
    urls = {}
    # Supplementary file URLs can be in gse.metadata["supplementary_file"]
    supp = gse.metadata.get("supplementary_file", [])
    for u in supp:
        # keep only the two expected files if present; otherwise download all
        for fn in FILES:
            if u.endswith("/" + fn) or u.endswith(fn):
                urls[fn] = u.replace("ftp://", "https://")
    # Fallback: if GEOparse didn't provide them, use direct
    if not urls:
        urls = direct_urls()
    return urls

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="data/raw", help="Output directory for downloaded files")
    ap.add_argument("--method", choices=["direct", "geoparse"], default="direct")
    ap.add_argument("--skip-existing", action="store_true", help="Do not re-download existing files")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    urls = direct_urls() if args.method == "direct" else geoparse_urls()

    print(f"[info] Downloading supplementary files for {GSE} using method={args.method}")
    for fn, url in urls.items():
        dest = outdir / fn
        if args.skip_existing and dest.exists():
            print(f"[skip] {dest} exists")
            continue
        print(f"[get] {url}")
        download_with_resume(url, dest)

    print("[done] downloads complete")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
