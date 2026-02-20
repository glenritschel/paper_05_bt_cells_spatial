import os
import zipfile
from pathlib import Path
from PIL import Image

core_dir = Path("reports/figures")
submission_dir = Path("submission_figures")
submission_dir.mkdir(exist_ok=True)

converted = []

# Convert all PNGs in reports/figures to 300 DPI JPG
for png in core_dir.rglob("*.png"):
    img = Image.open(png).convert("RGB")
    out = submission_dir / (png.stem + ".jpg")
    img.save(out, format="JPEG", dpi=(300,300), quality=95)
    converted.append(out)
    print("Converted:", out)

zip_path = Path("paper5_figures.zip")

with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
    for f in converted:
        zf.write(f, f.name)

print("Zip size:", os.path.getsize(zip_path), "bytes")

