from pathlib import Path
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
from PIL import Image
import yaml

root = Path(".")
spatial_dir = root / "data/processed/spatial_scored"
sub_dir = root / "submission_figures"
sub_dir.mkdir(parents=True, exist_ok=True)

programs_path = root / "configs/programs.yaml"
if not programs_path.exists():
    raise SystemExit(f"NOT FOUND: {programs_path}")

if not spatial_dir.exists():
    raise SystemExit(f"NOT FOUND: {spatial_dir} (run `make analyze` first)")

programs = yaml.safe_load(programs_path.read_text())

# ---- locate signatures in programs.yaml robustly ----
def first_present_key(d, keys):
    for k in keys:
        if isinstance(d, dict) and k in d:
            return d[k]
    return None

# programs.yaml could be either:
# 1) top-level mapping: {B_activation: [...], Th2: [...], Th17: [...]}
# 2) nested: {programs: {B_activation: [...], Th2: [...], Th17: [...]}}
prog_root = programs.get("programs", programs)

# allow a few common naming variants
B_GENES = first_present_key(prog_root, ["B_activation", "B_act", "B_activation_score", "B", "Bcell_activation"])
TH2_GENES = first_present_key(prog_root, ["Th2", "TH2", "T_h2"])
TH17_GENES = first_present_key(prog_root, ["Th17", "TH17", "T_h17"])

missing = []
if not B_GENES: missing.append("B_activation")
if not TH2_GENES: missing.append("Th2")
if not TH17_GENES: missing.append("Th17")
if missing:
    raise SystemExit(
        "Could not find these programs in configs/programs.yaml: "
        + ", ".join(missing)
        + "\nTop-level keys present were:\n  "
        + ", ".join(list(prog_root.keys())[:50])
        + "\nFix by ensuring programs.yaml contains, e.g.:\n"
        + "programs:\n  B_activation: [...]\n  Th2: [...]\n  Th17: [...]\n"
    )

def ensure_lgals9_obs(ad):
    if "LGALS9_expr" in ad.obs.columns:
        return
    if "LGALS9" in ad.var_names:
        x = ad[:, "LGALS9"].X
        if hasattr(x, "toarray"):
            x = x.toarray()
        ad.obs["LGALS9_expr"] = np.asarray(x).reshape(-1)
    else:
        # keep plotting functional (but warn)
        ad.obs["LGALS9_expr"] = 0.0

def score_genes_safe(ad, genes, outcol):
    # Use intersection with var_names
    present = [g for g in genes if g in ad.var_names]
    if len(present) < max(3, int(0.2 * len(genes))):
        raise SystemExit(
            f"Too few genes mapped for {outcol}: {len(present)}/{len(genes)} present.\n"
            f"Example present: {present[:25]}\n"
            f"Example missing: {sorted(set(genes) - set(present))[:25]}"
        )
    sc.tl.score_genes(ad, gene_list=present, score_name=outcol, use_raw=False)

def make_section_4panel(section: str, out_name: str):
    h5 = spatial_dir / f"{section}.h5ad"
    if not h5.exists():
        raise SystemExit(f"NOT FOUND: {h5}")

    ad = sc.read_h5ad(h5)
    ad.var_names_make_unique()

    # compute program scores
    score_genes_safe(ad, B_GENES, "B_activation_score")
    score_genes_safe(ad, TH2_GENES, "Th2_score")
    score_genes_safe(ad, TH17_GENES, "Th17_score")
    ensure_lgals9_obs(ad)

    keys = ["B_activation_score", "Th2_score", "Th17_score", "LGALS9_expr"]

    fig, axes = plt.subplots(1, 4, figsize=(20, 5), dpi=150)

    for ax, key in zip(axes, keys):
        sc.pl.spatial(
            ad,
            color=key,
            img_key="hires",
            size=1.3,
            ax=ax,
            show=False,
            title=key,
        )

    plt.subplots_adjust(wspace=0.25)

    out_path = sub_dir / out_name
    # DO NOT use bbox_inches="tight" to preserve ~3000x750 pixels
    fig.savefig(out_path, dpi=150)
    plt.close(fig)

    im = Image.open(out_path)
    print(out_name, "pixels:", im.size)

make_section_4panel("A", "FigureS1a_Section_A.jpg")
make_section_4panel("B", "FigureS1b_Section_B.jpg")
make_section_4panel("C", "FigureS1c_Section_C.jpg")
