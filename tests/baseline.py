"""Helpers for production-data baseline tests."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


BASELINE_PATH = Path(__file__).parent / "baselines" / "production_data.json"


def load_production_baseline() -> dict[str, Any]:
    with BASELINE_PATH.open(encoding="utf-8") as handle:
        return json.load(handle)
