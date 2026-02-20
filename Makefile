SHELL := /usr/bin/env bash
PY := python
CFG := configs/config.yaml

.PHONY: help setup download extract load notebook all clean

help:
	@echo "Targets:"
	@echo "  setup     - sanity check environment + create kernel"
	@echo "  download  - download GEO supplementary files"
	@echo "  extract   - extract RAW tar and verify Visium structure"
	@echo "  load      - load Visium sections into AnnData (.h5ad)"
	@echo "  notebook  - start Jupyter"
	@echo "  all       - download + extract + load"
	@echo "  clean     - remove interim/processed (keeps raw downloads)"

setup:
	@echo "[info] Python: $$($(PY) --version)"
	@$(PY) -c "import scanpy, squidpy; print('scanpy', scanpy.__version__, 'squidpy', squidpy.__version__)"
	@$(PY) -m ipykernel install --user --name paper05-bt-spatial --display-name "paper05-bt-spatial" >/dev/null 2>&1 || true
	@echo "[done] kernel installed (or already present)"

download:
	$(PY) scripts/00_download_geo.py --config $(CFG) --skip-existing

extract:
	$(PY) scripts/01_extract_and_verify.py --config $(CFG) $(if $(FORCE),--force,)

load:
	$(PY) scripts/02_load_visium.py --config $(CFG)

all: download extract load

notebook:
	jupyter notebook notebooks

preprocess:
	$(PY) scripts/03_preprocess_visium.py --config $(CFG)

score:
	$(PY) scripts/04_score_bt_markers.py --config $(CFG)

metrics:
	$(PY) scripts/05_bt_colocalization_metrics.py --config $(CFG)

analyze: preprocess score metrics

programs:
	$(PY) scripts/06_score_programs.py

corr:
	$(PY) scripts/07_bt_th2_th17_correlations.py

lgals9:
	$(PY) scripts/08_lgals9_distance_gradient.py

corr-immune:
	$(PY) scripts/07_bt_th2_th17_correlations.py --immune-only --immune-q 0.9

sens:
	$(PY) scripts/07_bt_th2_th17_correlations.py --immune-only --immune-q 0.85 --outcsv reports/b_th2_th17_correlations_q085.csv
	$(PY) scripts/07_bt_th2_th17_correlations.py --immune-only --immune-q 0.95 --outcsv reports/b_th2_th17_correlations_q095.csv

lgals9-immune:
	$(PY) scripts/08_lgals9_distance_gradient.py --immune-only --immune-q 0.9 --adjust --outcsv reports/lgals9_distance_effects_immune_q090.csv

composite:
	$(PY) scripts/10_make_composite_figure.py --section D --corr-csv reports/b_th2_th17_correlations.csv

clean:
	rm -rf data/interim/* data/processed/*
	rm -f data/interim/.outer_extract.stamp.json
	find data/interim -name ".extract.stamp.json" -delete 2>/dev/null || true
	@touch data/interim/.gitkeep data/processed/.gitkeep
	@echo "[done] cleaned interim/processed + removed extract stamps"

clean-processed:
	rm -rf data/processed/*
	@touch data/processed/.gitkeep
	@echo "[done] cleaned processed only"

