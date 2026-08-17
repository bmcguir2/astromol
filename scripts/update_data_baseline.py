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
    detection_rate_by_atoms_data,
    du_by_source_type_data,
    du_histogram_data,
    individual_source_data,
    kappa_histogram_data,
    mass_by_source_type_data,
    molecule_type_by_source_type_data,
    molecule_type_data,
    plot_du_bar_chart,
    plot_du_histogram,
    plot_molecule_type_by_source_enrichment_matrix,
    relative_du_by_source_type_data,
    rolling_rate_by_atoms_heatmap_data,
    scopes_by_year_data,
    source_type_data,
    wavelength_by_source_type_data,
)
from astromol.latex import (  # noqa: E402
    balanced_ism_table_columns,
    exgal_table_detections,
    exoplanet_table_detections,
    facility_table_entries,
    ice_table_detections,
    ppd_table_detections,
    rate_by_atoms_fits,
    rate_by_atoms_table_fragments,
    source_table_entries,
)
from astromol.slides import (  # noqa: E402
    build_molecule_slide_layout,
    build_ppd_detection_slide_layout,
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


def _axis_texts(figure: object, axes: object) -> list[str]:
    try:
        return [text.get_text() for text in axes.texts if text.get_text()]
    finally:
        import matplotlib.pyplot as plt

        plt.close(figure)


def _enrichment_matrix_text_prefix(data: object, count: int = 10) -> list[str]:
    figure, axes = plot_molecule_type_by_source_enrichment_matrix(data)
    return _axis_texts(figure, axes)[:count]


def build_regression_counts(
    db: Database,
    *,
    include_output_regressions: bool = True,
) -> dict[str, object]:
    """Return generated count expectations used by curation-sensitive tests."""
    view_2026 = CensusView.for_census(db, "2026")
    exgal_detections = exgal_table_detections(view_2026)
    exoplanet_detections = exoplanet_table_detections(view_2026)
    expanded_exoplanet_detections = exoplanet_table_detections(
        view_2026,
        include_isotopologues=True,
    )
    ice_detections = ice_table_detections(view_2026)
    secure_only_ice_detections = ice_table_detections(
        view_2026,
        include_tentative=False,
    )
    ppd_detections = ppd_table_detections(view_2026)
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
        "latex_exgal_table_2026": {
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
        "latex_ice_table_2026": {
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
        "latex_ppd_table_2026": {
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
    if not include_output_regressions or len(view_2026.ism_molecules()) < 10:
        return counts

    cumulative_detections = cumulative_detection_data(view_2026)
    cumulative_by_atoms = cumulative_by_atoms_data(view_2026)
    detection_rate_by_atoms = detection_rate_by_atoms_data(view_2026)
    rate_fits = rate_by_atoms_fits(view_2026)
    rate_fragments = rate_by_atoms_table_fragments(view_2026)
    du_histogram = du_histogram_data(view_2026)
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
    ppd_slide = build_ppd_detection_slide_layout(view_2026)

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
        "figures_detection_rate_by_atoms_2026": {
            "points": _detection_rate_points(detection_rate_by_atoms),
        },
        "figures_du_histogram_2026": {
            "max_du": float(du_histogram.max_du),
            "max_du_labels": _du_max_labels(du_histogram),
            "histogram_text": _axis_texts(
                *plot_du_histogram(du_histogram),
            ),
            "bar_text": _axis_texts(
                *plot_du_bar_chart(du_histogram),
            ),
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
            "mass_range": [
                round(value, 3) for value in mass_by_source.mass_range
            ],
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
            "enrichment_matrix_text_prefix": _enrichment_matrix_text_prefix(
                molecule_type_by_source,
            ),
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
        "latex_rate_by_atoms_table_2026": {
            "fit_values": _rate_fit_values(rate_fits),
            "table_rows": [
                line
                for line in rate_fragments["rates_by_atoms_table.tex"].splitlines()
                if "\t&\t" in line
            ],
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
            "balanced_molecule_font_pt": balanced_slide.molecule_font_pt,
            "balanced_group_molecule_counts": [
                len(group.molecules) for group in balanced_slide.groups
            ],
            "balanced_group_column_counts": [
                len(group.columns) for group in balanced_slide.groups
            ],
        },
        "slides_ppd_detection_slide_2026": {
            "total": int(ppd_slide.total),
            "count_text": f"{ppd_slide.total} Molecules",
            "molecule_font_pt": ppd_slide.molecule_font_pt,
            "group_labels": [group.spec.label for group in ppd_slide.groups],
            "group_molecule_counts": [
                len(group.molecules) for group in ppd_slide.groups
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
