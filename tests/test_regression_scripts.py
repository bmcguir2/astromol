from __future__ import annotations

from pathlib import Path
import runpy

import pytest


REGRESSION_SCRIPT_DIR = Path(__file__).parent / "regression_scripts"
REGRESSION_SCRIPTS = tuple(sorted(REGRESSION_SCRIPT_DIR.glob("*.py")))


@pytest.mark.parametrize(
    "script_path",
    REGRESSION_SCRIPTS,
    ids=lambda path: path.stem,
)
def test_regression_script(script_path: Path) -> None:
    runpy.run_path(str(script_path), run_name="__main__")
