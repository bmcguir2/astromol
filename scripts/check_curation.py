"""Run local verification checks for production curation changes."""

from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]

CHECKS = (
    (
        "Validate production data",
        [sys.executable, "-m", "astromol.validation"],
    ),
    (
        "Run curation-sensitive tests",
        [
            sys.executable,
            "-m",
            "pytest",
            "tests/test_load.py",
            "tests/test_validation.py",
            "tests/test_regression_scripts.py",
        ],
    ),
)


def _check_environment() -> dict[str, str]:
    env = os.environ.copy()
    env.setdefault("MPLBACKEND", "Agg")
    if "MPLCONFIGDIR" not in env:
        mpl_config = Path(tempfile.gettempdir()) / "astromol-mplconfig"
        mpl_config.mkdir(parents=True, exist_ok=True)
        env["MPLCONFIGDIR"] = str(mpl_config)
    return env


def main() -> int:
    env = _check_environment()
    for label, command in CHECKS:
        print(f"\n==> {label}", flush=True)
        print(" ".join(command), flush=True)
        completed = subprocess.run(command, cwd=ROOT, env=env)
        if completed.returncode:
            return completed.returncode
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
