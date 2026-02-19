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
	$(PY) scripts/01_extract_and_verify.py --config $(CFG)

load:
	$(PY) scripts/02_load_visium.py --config $(CFG)

all: download extract load

notebook:
	jupyter notebook notebooks

clean:
	rm -rf data/interim/* data/processed/*
	@touch data/interim/.gitkeep data/processed/.gitkeep
	@echo "[done] cleaned interim/processed"
