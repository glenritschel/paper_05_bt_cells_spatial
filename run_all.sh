#!/usr/bin/env bash
set -euo pipefail

python scripts/00_download_geo.py --outdir data/raw --skip-existing
python scripts/01_extract_and_verify.py --raw-tar data/raw/GSE249279_RAW.tar --outdir data/raw/extracted
python scripts/02_load_visium.py --indir data/raw/extracted --outdir data/processed/spatial

echo
echo "Next: jupyter notebook notebooks/01_qc_visium.ipynb"
