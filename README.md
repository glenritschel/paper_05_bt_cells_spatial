# GSE249279 (Ma et al., Nat Comm 2024) — Spatial Transcriptomics Setup Pipeline

This repo provides a **reproducible Python 3.11 pipeline** to download, unpack, verify, and load the GEO dataset **GSE249279** into **Scanpy/Squidpy**.

Dataset: *Systems-based identification of the Hippo pathway for promoting fibrotic mesenchymal differentiation in systemic sclerosis* (Ma et al., 2024).  
GEO Series page: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE249279

## What this pipeline does

1. Downloads GEO supplementary files:
   - `GSE249279_RAW.tar`
   - `GSE249279_scRNA-seq_merged_raw_counts.h5`
2. Extracts `GSE249279_RAW.tar` and **verifies** presence of Visium-style assets:
   - count matrix (`filtered_feature_bc_matrix.h5` or `matrix.mtx.gz` + `barcodes.tsv.gz` + `features.tsv.gz`)
   - `spatial/` folder (`tissue_hires_image.png`, `tissue_lowres_image.png`, `scalefactors_json.json`, `tissue_positions*`)
3. Loads each Visium sample into **AnnData** using `scanpy.read_visium()`
4. Saves processed `.h5ad` per Visium sample
5. Provides a notebook to confirm loading and generate basic QC plots

## Install

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python -m ipykernel install --user --name gse249279
```

## Run end-to-end

```bash
bash run_all.sh
```

Or step-by-step:

```bash
python scripts/00_download_geo.py --outdir data/raw
python scripts/01_extract_and_verify.py --raw-tar data/raw/GSE249279_RAW.tar --outdir data/raw/extracted
python scripts/02_load_visium.py --indir data/raw/extracted --outdir data/processed/spatial
```

Then open:

```bash
jupyter notebook notebooks/01_qc_visium.ipynb
```

## Outputs

- `data/raw/` — downloaded supplementary files
- `data/raw/extracted/` — extracted contents of `GSE249279_RAW.tar`
- `data/raw/manifest.json` — file inventory + detected Visium sample directories
- `data/processed/spatial/*.h5ad` — AnnData per tissue section
