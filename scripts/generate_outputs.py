"""Generate standard astromol output products from a source checkout."""

from __future__ import annotations

from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from astromol.outputs import main  # noqa: E402


if __name__ == "__main__":
    raise SystemExit(main())
