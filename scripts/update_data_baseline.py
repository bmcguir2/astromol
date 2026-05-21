"""Update the committed production-data test baseline."""

from __future__ import annotations

import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from astromol.database import Database  # noqa: E402
from astromol.census import CensusView  # noqa: E402
from astromol.latex import exoplanet_table_detections  # noqa: E402
from astromol.validation import validate_database  # noqa: E402


BASELINE_PATH = ROOT / "tests" / "baselines" / "production_data.json"


def build_regression_counts(db: Database) -> dict[str, object]:
    """Return generated count expectations used by curation-sensitive tests."""
    view_2026 = CensusView.for_census(db, "2026")
    exoplanet_detections = exoplanet_table_detections(view_2026)
    expanded_exoplanet_detections = exoplanet_table_detections(
        view_2026,
        include_isotopologues=True,
    )

    return {
        "census_view_2026": {
            "ism_molecules": len(view_2026.ism_molecules()),
            "ism_molecules_with_isotopologues": len(
                view_2026.ism_molecules(include_isotopologues=True)
            ),
            "exoplanet_molecules": len(view_2026.exoplanet_molecules()),
            "exoplanet_molecules_with_isotopologues": len(
                view_2026.exoplanet_molecules(include_isotopologues=True)
            ),
            "ppd_molecules": len(view_2026.ppd_molecules()),
            "ppd_molecules_with_isotopologues": len(
                view_2026.ppd_molecules(include_isotopologues=True)
            ),
        },
        "latex_exoplanet_table_2026": {
            "detections": len(exoplanet_detections),
            "linked_labels": len(
                {detection.molecule.label for detection in exoplanet_detections}
            ),
            "detections_with_isotopologues": len(expanded_exoplanet_detections),
            "isotopologue_detections": sum(
                detection.molecule.isotopologue_of is not None
                for detection in expanded_exoplanet_detections
            ),
        },
    }


def build_baseline(db: Database | None = None) -> dict[str, object]:
    db = db or Database()
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
        "schema_version": 2,
        "counts": {
            "telescopes": len(db.telescopes),
            "sources": len(db.sources),
            "molecules": len(db.molecules),
            "detections": len(db.detections),
        },
        "regression_counts": build_regression_counts(db),
        "validation_warnings": warnings,
    }


def write_baseline(baseline: dict[str, object], path: Path = BASELINE_PATH) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(baseline, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def main() -> int:
    baseline = build_baseline()
    write_baseline(baseline)
    print(f"Updated {BASELINE_PATH.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
