import os
import zipfile
from pathlib import Path
from PIL import Image

def open_img(p: Path) -> Image.Image:
    return Image.open(p).convert("RGB")

def save_jpeg(img: Image.Image, out: Path, dpi=300, quality=95):
    out.parent.mkdir(parents=True, exist_ok=True)
    img.save(out, format="JPEG", dpi=(dpi, dpi), quality=quality)

def grid(images, rows, cols, pad=20, bg=(255, 255, 255)):
    assert len(images) == rows * cols, f"Expected {rows*cols} images, got {len(images)}"
    ws = [im.size[0] for im in images]
    hs = [im.size[1] for im in images]
    w = max(ws)
    h = max(hs)
    W = cols * w + (cols + 1) * pad
    H = rows * h + (rows + 1) * pad
    canvas = Image.new("RGB", (W, H), bg)
    for i, im in enumerate(images):
        r = i // cols
        c = i % cols
        x = pad + c * (w + pad)
        y = pad + r * (h + pad)
        # center each image in its cell
        ox = x + (w - im.size[0]) // 2
        oy = y + (h - im.size[1]) // 2
        canvas.paste(im, (ox, oy))
    return canvas

root = Path(".")
fig_core = root / "reports/figures/core"
fig_comp = root / "reports/figures/composite"
out_dir = root / "submission_figures"
out_dir.mkdir(exist_ok=True)

# -------------------------
# Figure 1 (Section D 4-panel)
# Use the already-composed composite PNG
# -------------------------
f1_png = fig_comp / "D__4panel__Bact_Th2_Th17_LGALS9.png"
f1_jpg = out_dir / "Figure1_Spatial_Maps_SectionD.jpg"
if not f1_png.exists():
    raise SystemExit(f"NOT FOUND: {f1_png}")
save_jpeg(open_img(f1_png), f1_jpg)

# -------------------------
# Figure 2 (8-panel scatter plots)
# 4 sections x (Th2, Th17) -> 8 PNGs in reports/figures/core
# Arrange as 4 rows (A–D) x 2 cols (Th2, Th17)
# -------------------------
scatter_order = [
    fig_core / "A__Bact_vs_Th2__scatter.png",  fig_core / "A__Bact_vs_Th17__scatter.png",
    fig_core / "B__Bact_vs_Th2__scatter.png",  fig_core / "B__Bact_vs_Th17__scatter.png",
    fig_core / "C__Bact_vs_Th2__scatter.png",  fig_core / "C__Bact_vs_Th17__scatter.png",
    fig_core / "D__Bact_vs_Th2__scatter.png",  fig_core / "D__Bact_vs_Th17__scatter.png",
]
missing = [p for p in scatter_order if not p.exists()]
if missing:
    raise SystemExit("Missing scatter panels:\n" + "\n".join(str(p) for p in missing))

scatter_imgs = [open_img(p) for p in scatter_order]
f2_img = grid(scatter_imgs, rows=4, cols=2, pad=20)
f2_jpg = out_dir / "Figure2_Scatter_Plots.jpg"
save_jpeg(f2_img, f2_jpg)

# -------------------------
# Figure 3 (Summary barplot)
# Prefer the composite summary, which is already combined
# -------------------------
f3_png = fig_comp / "summary__adj_spearman__Bact_vs_Th2_Th17.png"
f3_jpg = out_dir / "Figure3_Summary_Barplot.jpg"
if not f3_png.exists():
    raise SystemExit(f"NOT FOUND: {f3_png}")
save_jpeg(open_img(f3_png), f3_jpg)

# -------------------------
# Figure S1 (Spatial maps A, B, C)
# Build from per-panel spatial PNGs in core:
# A/B/C each has 3 panels present? (You only listed D spatial panels in core)
# So we’ll build S1 from bt_scores + core if needed.
#
# Given your listing, only D has spatial program panels in core.
# Therefore we will use bt_scores as a surrogate if you intended those,
# OR (preferred) you should generate A/B/C spatial program panels.
#
# For now: build S1 using existing B/T score maps from bt_scores:
# For each section A/B/C: (B_score, T_score, Plasma_score, B_vs_T)
# If you want Bact/Th2/Th17/LGALS9 instead, tell me and I’ll give the exact generator.
# -------------------------
bt_scores = root / "reports/figures/bt_scores"
s1_panels = [
    bt_scores / "A__B_score.png", bt_scores / "A__T_score.png", bt_scores / "A__Plasma_score.png", bt_scores / "A__B_vs_T.png",
    bt_scores / "B__B_score.png", bt_scores / "B__T_score.png", bt_scores / "B__Plasma_score.png", bt_scores / "B__B_vs_T.png",
    bt_scores / "C__B_score.png", bt_scores / "C__T_score.png", bt_scores / "C__Plasma_score.png", bt_scores / "C__B_vs_T.png",
]
missing = [p for p in s1_panels if not p.exists()]
if missing:
    raise SystemExit("Missing S1 panels:\n" + "\n".join(str(p) for p in missing))

s1_imgs = [open_img(p) for p in s1_panels]
fS1_img = grid(s1_imgs, rows=3, cols=4, pad=20)
fS1_jpg = out_dir / "FigureS1_Spatial_Maps_ABC.jpg"
save_jpeg(fS1_img, fS1_jpg)

# -------------------------
# Figure S2 (Sensitivity analysis)
# You currently do NOT list q=0.85/q=0.95 scatter images in reports/figures.
# So we will create S2 as a 1x2 grid of the two sensitivity CSV barplots IF present,
# but those aren’t listed either.
#
# Minimal: reuse the same scatter layout but label “q085” and “q095” requires those images.
# Since they are not present, we will create S2 from existing correlation CSVs would need plotting.
# To keep this pure “convert existing PNGs” approach, we will instead package placeholders as NOT FOUND.
# -------------------------
# If you DO have these files, set them here:
s2_left = root / "reports/figures/sensitivity/q085__scatter_grid.png"
s2_right = root / "reports/figures/sensitivity/q095__scatter_grid.png"

if s2_left.exists() and s2_right.exists():
    s2_imgs = [open_img(s2_left), open_img(s2_right)]
    fS2_img = grid(s2_imgs, rows=1, cols=2, pad=20)
    fS2_jpg = out_dir / "FigureS2_Sensitivity_Analysis.jpg"
    save_jpeg(fS2_img, fS2_jpg)
else:
    print("NOTE: Sensitivity scatter grids not found; skipping FigureS2_Sensitivity_Analysis.jpg")
    print("Expected (if generated):")
    print(" -", s2_left)
    print(" -", s2_right)

# -------------------------
# ZIP
# -------------------------
zip_path = root / "paper5_figures.zip"

figure_files = [
    f1_jpg,
    f2_jpg,
    f3_jpg,
    fS1_jpg,
]

# include S2 only if created
s2_jpg = out_dir / "FigureS2_Sensitivity_Analysis.jpg"
if s2_jpg.exists():
    figure_files.append(s2_jpg)

with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
    for f in figure_files:
        if f.exists():
            zf.write(f, f.name)
            print("Added:", f)
        else:
            print("NOT FOUND:", f)

print("Zip size:", os.path.getsize(zip_path), "bytes")
print("Wrote:", zip_path)

