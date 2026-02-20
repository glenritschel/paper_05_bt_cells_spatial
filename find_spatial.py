from pathlib import Path
import scanpy as sc
import numpy as np

root = Path(".")
spatial_dir = root / "data/processed/spatial_scored"
if not spatial_dir.exists():
    raise SystemExit(f"Expected scored h5ads at {spatial_dir} but folder not found.")

ad = sc.read_h5ad(spatial_dir / "A.h5ad")
print("Loaded:", spatial_dir / "A.h5ad")
print("obs columns (first 50):", list(ad.obs.columns)[:50])
print("score cols present:", [c for c in ad.obs.columns if "score" in c.lower()][:50])
