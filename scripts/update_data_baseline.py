"""Update the committed production-data test baseline."""

from __future__ import annotations

import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from astromol.census import CensusView  # noqa: E402
from astromol.database import Database  # noqa: E402
from astromol.figures import (  # noqa: E402
    DU_BY_SOURCE_TYPE_BOXPLOT_ORDER,
    cumulative_by_atoms_data,
    cumulative_detection_data,
    du_by_source_type_data,
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
    balanced_ism_table_columns,
    exoplanet_table_detections,
    facility_table_entries,
    source_table_entries,
)
from astromol.slides import build_molecule_slide_layout  # noqa: E402
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


def _ring_categories_by_count(data: object) -> tuple[object, ...]:
    return tuple(
        category
        for _, category in sorted(
            enumerate(data.categories),
            key=lambda indexed_category: (
                -indexed_category[1].count,
                indexed_category[0],
            ),
        )
    )


def _ring_text(data: object, *, label_attr: str = "label") -> list[str]:
    text = []
    for category in _ring_categories_by_count(data):
        text.extend([getattr(category, label_attr), f"{category.percent:.1f}%"])
    return text


def _type_ring_text(data: object) -> list[str]:
    categories = _ring_categories_by_count(data)
    return [
        category.plural_label
        for category in categories
    ] + [
        f"{category.percent:.1f}%"
        for category in categories
    ]


def _source_boxplot_n_labels(
    data: object,
    value_attr: str,
    *,
    include_negative_values: bool = True,
) -> list[str]:
    categories_by_key = {category.key: category for category in data.categories}
    labels = []
    for key in DU_BY_SOURCE_TYPE_BOXPLOT_ORDER:
        category = categories_by_key.get(key)
        if category is None:
            continue
        values = list(getattr(category, value_attr))
        if not include_negative_values:
            values = [value for value in values if value >= 0]
        if values:
            labels.append(f"n={len(values)}")
    return labels


def _wavelength_stacked_n_labels(data: object) -> list[str]:
    categories_by_key = {category.key: category for category in data.categories}
    labels = []
    for key in DU_BY_SOURCE_TYPE_BOXPLOT_ORDER:
        category = categories_by_key.get(key)
        if category is None:
            continue
        total = sum(
            category.count_for_display_wavelength(wavelength)
            for wavelength in ("cm", "mm", "sub-mm", "IR", "UV-Vis")
        )
        labels.append(f"n={total}")
    return labels


def build_regression_counts(
    db: Database,
    *,
    include_output_regressions: bool = True,
) -> dict[str, object]:
    """Return generated count expectations used by curation-sensitive tests."""
    view_2026 = CensusView.for_census(db, "2026")
    exoplanet_detections = exoplanet_table_detections(view_2026)
    expanded_exoplanet_detections = exoplanet_table_detections(
        view_2026,
        include_isotopologues=True,
    )
    counts = {
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
    if not include_output_regressions or len(view_2026.ism_molecules()) < 10:
        return counts

    cumulative_detections = cumulative_detection_data(view_2026)
    cumulative_by_atoms = cumulative_by_atoms_data(view_2026)
    source_type = source_type_data(view_2026)
    molecule_type = molecule_type_data(view_2026)
    individual_source = individual_source_data(view_2026)
    mass_by_source = mass_by_source_type_data(view_2026)
    du_by_source = du_by_source_type_data(view_2026)
    relative_du_by_source = relative_du_by_source_type_data(view_2026)
    kappa_histogram = kappa_histogram_data(view_2026)
    molecule_type_by_source = molecule_type_by_source_type_data(view_2026)
    wavelength_by_source = wavelength_by_source_type_data(view_2026)
    scopes_by_year = scopes_by_year_data(view_2026)
    yebes_row = _scope_rows_by_label(scopes_by_year)["Yebes 40-m"]
    facility_entries = facility_table_entries(view_2026)
    source_entries = source_table_entries(view_2026)
    balanced_ism_groups = balanced_ism_table_columns(view_2026)
    legacy_slide = build_molecule_slide_layout(view_2026)
    balanced_slide = build_molecule_slide_layout(view_2026, profile="balanced")

    counts.update({
        "figures_cumulative_detections_2026": {
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
        "figures_cumulative_by_atoms_2026": {
            "total": int(cumulative_by_atoms.total),
            "final_counts": _final_counts_by_label(cumulative_by_atoms),
            "rolling_rate_last_column": [
                float(value)
                for value in rolling_rate_by_atoms_heatmap_data(view_2026).matrix[:, -1]
            ],
        },
        "figures_source_type_2026": {
            "molecule_count": int(source_type.molecule_count),
            "counts": source_type.counts,
            "ring_text": _ring_text(source_type),
        },
        "figures_molecule_type_2026": {
            "molecule_count": int(molecule_type.molecule_count),
            "counts": molecule_type.counts,
            "ring_text": _type_ring_text(molecule_type),
        },
        "figures_individual_source_2026": {
            "molecule_count": int(individual_source.molecule_count),
            "counts": individual_source.counts,
            "ring_text": _ring_text(individual_source),
        },
        "figures_mass_by_source_type_2026": {
            "molecule_count": int(mass_by_source.molecule_count),
            "counts": mass_by_source.counts,
            "boxplot_n_labels": _source_boxplot_n_labels(
                mass_by_source,
                "masses",
            ),
        },
        "figures_du_by_source_type_2026": {
            "molecule_count": int(du_by_source.molecule_count),
            "counts": du_by_source.counts,
            "boxplot_n_labels": _source_boxplot_n_labels(
                du_by_source,
                "values",
                include_negative_values=False,
            ),
        },
        "figures_relative_du_by_source_type_2026": {
            "molecule_count": int(relative_du_by_source.molecule_count),
            "counts": relative_du_by_source.counts,
            "boxplot_n_labels": _source_boxplot_n_labels(
                relative_du_by_source,
                "values",
                include_negative_values=False,
            ),
        },
        "figures_kappas_2026": {
            "molecule_count": int(kappa_histogram.molecule_count),
            "histogram_max": int(kappa_histogram.histogram_counts().max()),
        },
        "figures_molecule_type_by_source_type_2026": {
            "molecule_count": int(molecule_type_by_source.molecule_count),
            "counts": molecule_type_by_source.counts,
            "source_counts": molecule_type_by_source.source_counts,
            "overall_type_counts": molecule_type_by_source.overall_type_counts,
        },
        "figures_wavelength_by_source_type_2026": {
            "molecule_count": int(wavelength_by_source.molecule_count),
            "counts": wavelength_by_source.counts,
            "stacked_n_labels": _wavelength_stacked_n_labels(wavelength_by_source),
        },
        "figures_scopes_by_year_2026": {
            "yebes_final_count": int(yebes_row.final_count),
            "yebes_rate_rounded": round(yebes_row.rate, 1),
        },
        "latex_facility_table_2026": {
            "entries": len(facility_entries),
            "credited_detections": sum(count for _, count in facility_entries),
            "top_entries": facility_entries[:5],
        },
        "latex_source_table_2026": {
            "entries": len(source_entries),
            "credited_detections": sum(count for _, count in source_entries),
            "top_entries": source_entries[:5],
        },
        "latex_ism_tables_2026": {
            "balanced_group_column_counts": [
                len(group) for group in balanced_ism_groups
            ],
            "linked_labels": len(view_2026.ism_molecules()),
        },
        "slides_molecule_slide_2026": {
            "total": int(legacy_slide.total),
            "count_text": f"{legacy_slide.total} Molecules",
            "legacy_warning_count": len(legacy_slide.warnings),
            "balanced_group_molecule_counts": [
                len(group.molecules) for group in balanced_slide.groups
            ],
            "balanced_group_column_counts": [
                len(group.columns) for group in balanced_slide.groups
            ],
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
    baseline = build_baseline()
    write_baseline(baseline)
    print(f"Updated {BASELINE_PATH.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
