from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict

import yaml

def read_yaml(path: str | Path) -> Dict[str, Any]:
    p = Path(path)
    with p.open("r") as f:
        return yaml.safe_load(f)

def write_json(obj: Any, path: str | Path) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w") as f:
        json.dump(obj, f, indent=2)

def ensure_dir(path: str | Path) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p
