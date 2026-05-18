"""Update the committed production-data test baseline."""

from __future__ import annotations

import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from astromol.database import Database  # noqa: E402
from astromol.validation import validate_database  # noqa: E402


BASELINE_PATH = ROOT / "tests" / "baselines" / "production_data.json"


def build_baseline() -> dict[str, object]:
    db = Database()
    report = validate_database(db)

    warnings = [
        {
            "severity": issue.severity,
            "code": issue.code,
            "record": issue.record,
            "message": issue.message,
        }
        for issue in report.warnings
    ]
    warnings.sort(key=lambda issue: (issue["code"], issue["record"], issue["message"]))

    return {
        "schema_version": 1,
        "counts": {
            "references": len(db.refs),
            "telescopes": len(db.telescopes),
            "sources": len(db.sources),
            "molecules": len(db.molecules),
            "detections": len(db.detections),
        },
        "validation_warnings": warnings,
    }


def main() -> int:
    baseline = build_baseline()
    BASELINE_PATH.parent.mkdir(parents=True, exist_ok=True)
    BASELINE_PATH.write_text(
        json.dumps(baseline, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(f"Updated {BASELINE_PATH.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
