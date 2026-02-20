#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import yaml
import scanpy as sc

from src.paper05.utils import ensure_dir


def read_programs(p: Path) -> dict[str, list[str]]:
    d = yaml.safe_load(p.read_text())
    return d["programs"]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", default="data/processed/spatial_prepped")
    ap.add_argument("--outdir", default="data/processed/spatial_programs")
    ap.add_argument("--programs", default="configs/programs.yaml")
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()

    programs = read_programs(Path(args.programs))
    indir = Path(args.indir)
    outdir = ensure_dir(args.outdir)

    paths = sorted(indir.glob("*.h5ad"))
    if not paths:
        raise SystemExit(f"No inputs found: {indir}/*.h5ad")

    for p in paths:
        out = outdir / p.name
        if out.exists() and not args.overwrite:
            print(f"[skip] {out}")
            continue

        print(f"[score] {p.name} -> {out.name}")
        ad = sc.read_h5ad(p)
        ad.var_names_make_unique()

        present = set(ad.var_names)
        ad.uns.setdefault("program_scoring", {})

        for prog, genes in programs.items():
            genes_present = [g for g in genes if g in present]
            ad.uns["program_scoring"][prog] = {
                "requested": genes,
                "present": genes_present,
                "n_present": len(genes_present),
            }
            if len(genes_present) < 3:
                print(f"  [warn] {p.stem} {prog}: only {len(genes_present)} genes present")
            sc.tl.score_genes(ad, genes_present, score_name=f"{prog}_score", use_raw=False)

        ad.write_h5ad(out)

    print(f"[done] wrote {len(paths)} files to {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

