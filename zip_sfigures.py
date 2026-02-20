import os, zipfile
from pathlib import Path

figure_files = [
  "submission_figures/Figure1_Spatial_Maps_SectionD.jpg",
  "submission_figures/Figure2_Scatter_Plots.jpg",
  "submission_figures/Figure3_Summary_Barplot.jpg",
  "submission_figures/FigureS1a_Section_A.jpg",
  "submission_figures/FigureS1b_Section_B.jpg",
  "submission_figures/FigureS1c_Section_C.jpg",
]

out = Path("paper5_figures.zip")
with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as zf:
    for f in figure_files:
        if os.path.exists(f):
            zf.write(f, Path(f).name)
            print("Added:", f)
        else:
            print("NOT FOUND:", f)

print("Zip size:", out.stat().st_size, "bytes")
