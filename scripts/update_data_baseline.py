"""Update the committed production-data test baseline."""

from __future__ import annotations

import json
import argparse
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from astromol.census import CensusView  # noqa: E402
from astromol.database import Database  # noqa: E402
from astromol.figures import (  # noqa: E402
    cumulative_by_atoms_data,
    cumulative_detection_data,
    detection_rate_by_atoms_data,
    du_by_source_type_data,
    du_histogram_data,
    individual_source_data,
    kappa_histogram_data,
    mass_by_source_type_data,
    molecule_type_by_source_type_data,
    molecule_type_data,
    relative_du_by_source_type_data,
    rolling_rate_by_atoms_heatmap_data,
    scopes_by_year_data,
    source_type_data,
    wavelength_by_source_type_data,
)
from astromol.latex import (  # noqa: E402
    exgal_table_detections,
    exoplanet_table_detections,
    facility_table_entries,
    ice_table_detections,
    ppd_table_detections,
    rate_by_atoms_fits,
    source_table_entries,
)
from astromol.validation import validate_database  # noqa: E402


BASELINE_PATH = ROOT / "tests" / "baselines" / "production_data.json"


def _cumulative_count_at_year(data: object, year: int) -> int:
    for data_year, count in zip(data.years, data.counts):
        if int(data_year) == year:
            return int(count)
    raise ValueError(f"No cumulative count for year {year}")


def _final_counts_by_label(data: object) -> dict[str, int]:
    return {series.label: int(series.final_count) for series in data.series}


def _scope_rows_by_label(data: object) -> dict[str, object]:
    return {series.label: series for series in data.series}


def _detection_rate_points(data: object) -> dict[str, dict[str, float | int]]:
    return {
        point.label: {
            "count": int(point.count),
            "first_year": int(point.first_year),
            "rate": float(point.rate),
        }
        for point in data.points
    }


def _rate_fit_values(fits: object) -> dict[str, list[float | int]]:
    return {
        fit.label: [round(fit.slope, 2), round(fit.r_squared, 2), fit.onset_year]
        for fit in fits
    }


def _du_max_labels(data: object) -> list[str]:
    return [
        label
        for label, value in zip(data.molecule_labels, data.values)
        if value == data.max_du
    ]


def baseline_end_year(db: Database) -> int:
    """Stable endpoint for reviewed data snapshots, independent of rebuild date."""
    return max((d.year for d in db.detections if d.status == "secure"), default=2021)


def build_regression_counts(
    db: Database,
    *,
    include_output_regressions: bool = True,
) -> dict[str, object]:
    """Return generated count expectations used by curation-sensitive tests."""
    view_current = CensusView.current(db, end_year=baseline_end_year(db))
    exgal_detections = exgal_table_detections(view_current)
    exoplanet_detections = exoplanet_table_detections(view_current)
    expanded_exoplanet_detections = exoplanet_table_detections(
        view_current,
        include_isotopologues=True,
    )
    ice_detections = ice_table_detections(view_current)
    secure_only_ice_detections = ice_table_detections(
        view_current,
        include_tentative=False,
    )
    ppd_detections = ppd_table_detections(view_current)
    counts = {
        "census_view_current": {
            "ism_molecules": len(view_current.ism_molecules()),
            "ism_molecules_with_isotopologues": len(
                view_current.ism_molecules(include_isotopologues=True)
            ),
            "exoplanet_molecules": len(view_current.exoplanet_molecules()),
            "exoplanet_molecules_with_isotopologues": len(
                view_current.exoplanet_molecules(include_isotopologues=True)
            ),
            "ppd_molecules": len(view_current.ppd_molecules()),
            "ppd_molecules_with_isotopologues": len(
                view_current.ppd_molecules(include_isotopologues=True)
            ),
        },
        "latex_exgal_table_current": {
            "detections": len(exgal_detections),
            "linked_labels": len(
                {detection.molecule.label for detection in exgal_detections}
            ),
            "secure_detections": sum(
                detection.status == "secure" for detection in exgal_detections
            ),
            "tentative_detections": sum(
                detection.status == "tentative" for detection in exgal_detections
            ),
            "observation_references": len(
                {
                    ref.bibcode
                    for detection in exgal_detections
                    for ref in detection.refs.get("observation", [])
                }
            ),
        },
        "latex_exoplanet_table_current": {
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
        "latex_ice_table_current": {
            "detections": len(ice_detections),
            "linked_labels": len(
                {detection.molecule.label for detection in ice_detections}
            ),
            "secure_detections": sum(
                detection.status == "secure" for detection in ice_detections
            ),
            "tentative_detections": sum(
                detection.status == "tentative" for detection in ice_detections
            ),
            "secure_only_detections": len(secure_only_ice_detections),
            "observation_references": len(
                {
                    ref.bibcode
                    for detection in ice_detections
                    for ref in detection.refs.get("observation", [])
                }
            ),
        },
        "latex_ppd_table_current": {
            "detections": len(ppd_detections),
            "linked_labels": len(
                {detection.molecule.label for detection in ppd_detections}
            ),
            "isotopologue_detections": sum(
                detection.molecule.isotopologue_of is not None
                for detection in ppd_detections
            ),
        },
    }
    if not include_output_regressions or len(view_current.ism_molecules()) < 10:
        return counts

    cumulative_detections = cumulative_detection_data(view_current)
    cumulative_by_atoms = cumulative_by_atoms_data(view_current)
    detection_rate_by_atoms = detection_rate_by_atoms_data(view_current)
    rate_fits = rate_by_atoms_fits(view_current)
    du_histogram = du_histogram_data(view_current)
    source_type = source_type_data(view_current)
    molecule_type = molecule_type_data(view_current)
    individual_source = individual_source_data(view_current)
    mass_by_source = mass_by_source_type_data(view_current)
    du_by_source = du_by_source_type_data(view_current)
    relative_du_by_source = relative_du_by_source_type_data(view_current)
    kappa_histogram = kappa_histogram_data(view_current)
    molecule_type_by_source = molecule_type_by_source_type_data(view_current)
    wavelength_by_source = wavelength_by_source_type_data(view_current)
    scopes_by_year = scopes_by_year_data(view_current)
    yebes_row = _scope_rows_by_label(scopes_by_year)["Yebes 40-m"]
    facility_entries = facility_table_entries(view_current)
    source_entries = source_table_entries(view_current)

    counts.update({
        "figures_cumulative_detections_current": {
            "total": int(cumulative_detections.total),
            "first_detection_years": len(cumulative_detections.first_detection_years),
            "counts_by_year": {
                str(year): _cumulative_count_at_year(cumulative_detections, year)
                for year in (2019, 2020, 2021, 2024, 2026)
            },
            "trend_slopes": {
                trend.label: trend.slope
                for trend in cumulative_detections.trends
            },
        },
        "figures_cumulative_by_atoms_current": {
            "total": int(cumulative_by_atoms.total),
            "final_counts": _final_counts_by_label(cumulative_by_atoms),
            "rolling_rate_last_column": [
                float(value)
                for value in rolling_rate_by_atoms_heatmap_data(view_current).matrix[:, -1]
            ],
        },
        "figures_detection_rate_by_atoms_current": {
            "points": _detection_rate_points(detection_rate_by_atoms),
        },
        "figures_du_histogram_current": {
            "max_du": float(du_histogram.max_du),
            "max_du_labels": _du_max_labels(du_histogram),
        },
        "figures_source_type_current": {
            "molecule_count": int(source_type.molecule_count),
            "counts": source_type.counts,
        },
        "figures_molecule_type_current": {
            "molecule_count": int(molecule_type.molecule_count),
            "counts": molecule_type.counts,
        },
        "figures_individual_source_current": {
            "molecule_count": int(individual_source.molecule_count),
            "counts": individual_source.counts,
        },
        "figures_mass_by_source_type_current": {
            "molecule_count": int(mass_by_source.molecule_count),
            "counts": mass_by_source.counts,
            "mass_range": [
                round(value, 3) for value in mass_by_source.mass_range
            ],
        },
        "figures_du_by_source_type_current": {
            "molecule_count": int(du_by_source.molecule_count),
            "counts": du_by_source.counts,
        },
        "figures_relative_du_by_source_type_current": {
            "molecule_count": int(relative_du_by_source.molecule_count),
            "counts": relative_du_by_source.counts,
        },
        "figures_kappas_current": {
            "molecule_count": int(kappa_histogram.molecule_count),
            "histogram_max": int(kappa_histogram.histogram_counts().max()),
        },
        "figures_molecule_type_by_source_type_current": {
            "molecule_count": int(molecule_type_by_source.molecule_count),
            "counts": molecule_type_by_source.counts,
            "source_counts": molecule_type_by_source.source_counts,
            "overall_type_counts": molecule_type_by_source.overall_type_counts,
        },
        "figures_wavelength_by_source_type_current": {
            "molecule_count": int(wavelength_by_source.molecule_count),
            "counts": wavelength_by_source.counts,
        },
        "figures_scopes_by_year_current": {
            "yebes_final_count": int(yebes_row.final_count),
            "yebes_rate_rounded": round(yebes_row.rate, 1),
        },
        "latex_facility_table_current": {
            "entries": len(facility_entries),
            "credited_detections": sum(count for _, count in facility_entries),
            "top_entries": facility_entries[:5],
        },
        "latex_source_table_current": {
            "entries": len(source_entries),
            "credited_detections": sum(count for _, count in source_entries),
            "top_entries": source_entries[:5],
        },
        "latex_rate_by_atoms_table_current": {
            "fit_values": _rate_fit_values(rate_fits),
        },
        "latex_ism_tables_current": {
            "linked_labels": len(view_current.ism_molecules()),
        },
        "slides_molecule_slide_current": {
            "total": len(view_current.ism_molecules()),
        },
        "slides_ppd_detection_slide_current": {
            "total": len(view_current.ppd_molecules(include_isotopologues=True)),
        },
    })
    return counts


def build_baseline(
    db: Database | None = None,
    *,
    include_output_regressions: bool = True,
) -> dict[str, object]:
    db = db or Database()
    report = validate_database(db)
    report.raise_for_errors()

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
        "schema_version": 3,
        "analysis_end_year": baseline_end_year(db),
        "membership": {
            name: sorted(d.id for d in view.detections(include_tentative=True, include_disputed=True, include_isotopologues=True))
            for name, view in [("2018", CensusView.for_census(db, "2018")),
                               ("2021", CensusView.for_census(db, "2021")),
                               ("current", CensusView.current(db))]
        },
        "counts": {
            "telescopes": len(db.telescopes),
            "sources": len(db.sources),
            "molecules": len(db.molecules),
            "detections": len(db.detections),
        },
        "regression_counts": build_regression_counts(
            db,
            include_output_regressions=include_output_regressions,
        ),
        "validation_warnings": warnings,
    }


def write_baseline(baseline: dict[str, object], path: Path = BASELINE_PATH) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(baseline, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def main() -> int:
    argparse.ArgumentParser(description=__doc__).parse_args()
    baseline = build_baseline()
    write_baseline(baseline)
    print(f"Updated {BASELINE_PATH.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
