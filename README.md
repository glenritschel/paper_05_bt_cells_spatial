# Paper 05 — SSc Skin Spatial Transcriptomics (B↔T focus)

Reproducible **micromamba** pipeline to download, verify, and load **GSE249279** (Ma et al., Nat Comm 2024) Visium spatial data into **Scanpy/Squidpy**, and generate basic QC plots.

## Quickstart

```bash
micromamba env create -f env/environment.yml
micromamba activate paper05-bt-spatial
make setup
make all
make notebook
```

## Layout

- `scripts/` — pipeline steps (download, extract/verify, load)
- `src/paper05/` — small utilities
- `configs/` — config YAML (paths)
- `notebooks/` — QC notebook
- `data/` + `reports/` — outputs (gitignored; `.gitkeep` kept)
