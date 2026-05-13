"""Generate standard astromol output bundles.

This module powers the GitHub Actions workflow that publishes the latest
standard figures and PowerPoint slides. It intentionally uses the public table,
figure, and slide APIs so automated products exercise the same code paths that
users call from notebooks or local scripts.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
from html import escape
from pathlib import Path
import argparse
import os
import shutil
import zipfile

os.environ.setdefault("MPLBACKEND", "Agg")

from .census import CensusView
from .database import Database
from .figures import (
    cumulative_by_atoms_data,
    cumulative_detection_data,
    detection_rate_by_atoms_data,
    du_by_source_type_data,
    du_histogram_data,
    facility_share_data,
    individual_source_data,
    kappa_histogram_data,
    mass_by_source_type_data,
    mass_by_wavelength_data,
    molecule_type_by_source_type_data,
    molecule_type_data,
    molecules_by_wavelength_atoms_data,
    periodic_heatmap_data,
    relative_du_by_source_type_data,
    rolling_rate_by_atoms_heatmap_data,
    scopes_by_year_data,
    source_type_data,
    wavelength_by_source_type_data,
    write_cumulative_by_atoms_plot,
    write_cumulative_detections_plot,
    write_detection_rate_by_atoms_comparison_plot,
    write_detection_rate_by_atoms_plot,
    write_du_bar_chart,
    write_du_by_source_type_boxplot,
    write_facility_share_bars_plot,
    write_individual_source_pie_chart,
    write_kappa_histogram,
    write_mass_by_source_type_boxplot,
    write_mass_by_wavelength_boxplot,
    write_mass_by_wavelength_plot,
    write_molecule_type_by_source_enrichment_matrix,
    write_molecules_by_wavelength_atoms_bubble_heatmap,
    write_periodic_heatmap,
    write_relative_du_by_source_type_boxplot,
    write_rolling_rate_by_atoms_heatmap,
    write_scopes_by_year_plot,
    write_source_pie_chart,
    write_stacked_cumulative_by_atoms_plot,
    write_type_pie_chart,
    write_wavelength_by_source_type_stacked_bar,
)
from .latex import (
    write_exgal_table,
    write_exoplanet_table,
    write_facility_table,
    write_ice_table,
    write_ism_tables,
    write_ppd_table,
    write_rate_by_atoms_table,
    write_scalar_fragments,
    write_source_table,
)
from .slides import (
    build_molecule_slide_layout,
    build_ppd_detection_slide_layout,
    write_molecule_slide,
    write_molecule_slide_report,
    write_ppd_detection_slide,
)


DEFAULT_OUTPUT_DIR = Path("build") / "astromol_outputs"
DEFAULT_FORMATS = ("png", "pdf")


@dataclass(frozen=True)
class GeneratedProduct:
    """One generated output file for the static product index."""

    label: str
    path: Path
    description: str


def view_from_choice(db: Database, choice: str) -> CensusView:
    """Return a census/current view from a command-line choice."""
    if choice == "current":
        return CensusView.current(db)
    return CensusView.for_census(db, choice)


def generate_standard_outputs(
    output_dir: str | Path = DEFAULT_OUTPUT_DIR,
    *,
    view_choice: str = "current",
    formats: tuple[str, ...] = DEFAULT_FORMATS,
    clean: bool = True,
    include_tables: bool = True,
    include_slides: bool = True,
) -> list[GeneratedProduct]:
    """Generate standard latest astromol products and return their paths."""
    output_dir = Path(output_dir)
    if clean and output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    db = Database()
    view = view_from_choice(db, view_choice)
    view_2021 = CensusView.for_census(db, "2021")
    products: list[GeneratedProduct] = []

    products.extend(
        _generate_figures(
            view,
            view_2021=view_2021,
            output_dir=output_dir,
            view_choice=view_choice,
            formats=formats,
        )
    )
    if include_slides:
        products.extend(
            _generate_slides(
                view,
                output_dir=output_dir,
                view_choice=view_choice,
            )
        )
    if include_tables:
        products.extend(
            _generate_tables(
                view,
                output_dir=output_dir,
            )
        )

    zip_path = _write_zip(output_dir)
    products.append(
        GeneratedProduct(
            label="Complete output bundle",
            path=zip_path,
            description="Zip archive containing all generated standard products.",
        )
    )
    _write_index(
        output_dir,
        products,
        view_choice=view_choice,
        db=db,
    )
    return products


def _generate_figures(
    view: CensusView,
    *,
    view_2021: CensusView,
    output_dir: Path,
    view_choice: str,
    formats: tuple[str, ...],
) -> list[GeneratedProduct]:
    """Generate the standard manuscript-facing figure set."""
    figure_dir = output_dir / "figures"
    products: list[GeneratedProduct] = []

    by_atoms = cumulative_by_atoms_data(view)
    du_source = du_by_source_type_data(view)
    relative_du_source = relative_du_by_source_type_data(view)
    mass_source = mass_by_source_type_data(view)
    wavelength_source = wavelength_by_source_type_data(view)
    rate_data = detection_rate_by_atoms_data(view)
    rate_2021 = detection_rate_by_atoms_data(view_2021)

    figure_specs = [
        (
            "cumulative_detections",
            "Cumulative ISM/CSM detections",
            write_cumulative_detections_plot,
            cumulative_detection_data(view),
            {},
        ),
        (
            "cumulative_by_atoms",
            "Cumulative detections by atom-count category",
            write_cumulative_by_atoms_plot,
            by_atoms,
            {},
        ),
        (
            "cumulative_by_atoms_stacked",
            "Stacked cumulative detections by atom-count category",
            write_stacked_cumulative_by_atoms_plot,
            by_atoms,
            {},
        ),
        (
            "rolling_rate_by_atoms",
            "Rolling detection-rate heatmap by atom-count category",
            write_rolling_rate_by_atoms_heatmap,
            rolling_rate_by_atoms_heatmap_data(view),
            {},
        ),
        (
            "periodic_heatmap",
            "Periodic-table heatmap of elements in detected molecules",
            write_periodic_heatmap,
            periodic_heatmap_data(view),
            {},
        ),
        (
            "du_bar_chart",
            "Degree-of-unsaturation exact-value bar chart",
            write_du_bar_chart,
            du_histogram_data(view),
            {},
        ),
        (
            "kappa_distribution",
            "Ray asymmetry-parameter distribution",
            write_kappa_histogram,
            kappa_histogram_data(view),
            {},
        ),
        (
            "molecule_type_pie",
            "Molecule-type pie chart",
            write_type_pie_chart,
            molecule_type_data(view),
            {},
        ),
        (
            "source_type_pie",
            "Source-type pie chart",
            write_source_pie_chart,
            source_type_data(view),
            {},
        ),
        (
            "individual_source_pie",
            "Individual-source contribution pie chart",
            write_individual_source_pie_chart,
            individual_source_data(view),
            {},
        ),
        (
            "molecule_type_source_enrichment",
            "Molecule-type enrichment by source type",
            write_molecule_type_by_source_enrichment_matrix,
            molecule_type_by_source_type_data(view),
            {},
        ),
        (
            "du_by_source_boxplot",
            "Degree of unsaturation by source type",
            write_du_by_source_type_boxplot,
            du_source,
            {},
        ),
        (
            "relative_du_by_source_boxplot",
            "Relative degree of unsaturation by source type",
            write_relative_du_by_source_type_boxplot,
            relative_du_source,
            {},
        ),
        (
            "mass_by_source_boxplot",
            "Molecular mass by source type",
            write_mass_by_source_type_boxplot,
            mass_source,
            {},
        ),
        (
            "wavelength_by_source_stacked_bar",
            "First-detection wavelength share by source type",
            write_wavelength_by_source_type_stacked_bar,
            wavelength_source,
            {},
        ),
        (
            "mass_by_wavelength_kde",
            "Molecular mass distribution by first-detection wavelength",
            write_mass_by_wavelength_plot,
            mass_by_wavelength_data(view, include_fullerenes=False),
            {},
        ),
        (
            "mass_by_wavelength_boxplot",
            "Molecular mass by first-detection wavelength",
            write_mass_by_wavelength_boxplot,
            mass_by_wavelength_data(view, include_fullerenes=False),
            {},
        ),
        (
            "wavelength_atoms_bubble_heatmap",
            "Detection wavelength by atom-count bubble heatmap",
            write_molecules_by_wavelength_atoms_bubble_heatmap,
            molecules_by_wavelength_atoms_data(view),
            {},
        ),
        (
            "detection_rate_by_atoms",
            "Detection rate by atom-count category",
            write_detection_rate_by_atoms_plot,
            rate_data,
            {},
        ),
        (
            "detection_rate_by_atoms_comparison",
            "Detection rate by atom-count category with 2021 baseline",
            write_detection_rate_by_atoms_comparison_plot,
            rate_data,
            {"baseline_data": rate_2021, "current_label": view_choice},
        ),
        (
            "facility_share_bars",
            "Facility contribution share by facility era",
            write_facility_share_bars_plot,
            facility_share_data(view),
            {},
        ),
        (
            "scopes_by_year",
            "Cumulative facility contributions by year",
            write_scopes_by_year_plot,
            scopes_by_year_data(view),
            {"style": "modern"},
        ),
    ]

    for name, label, writer, data, kwargs in figure_specs:
        for file_format in formats:
            path = figure_dir / file_format / f"{name}.{file_format}"
            path.parent.mkdir(parents=True, exist_ok=True)
            if "baseline_data" in kwargs:
                baseline_data = kwargs["baseline_data"]
                writer(
                    data,
                    baseline_data,
                    path,
                    current_label=kwargs.get("current_label", view_choice),
                )
            else:
                writer(data, path, **kwargs)
            products.append(
                GeneratedProduct(
                    label=f"{label} ({file_format.upper()})",
                    path=path,
                    description="Standard latest figure output.",
                )
            )
    return products


def _generate_slides(
    view: CensusView,
    *,
    output_dir: Path,
    view_choice: str,
) -> list[GeneratedProduct]:
    """Generate standard PowerPoint slides and layout reports."""
    slide_dir = output_dir / "slides"
    slide_dir.mkdir(parents=True, exist_ok=True)
    products: list[GeneratedProduct] = []

    ism_slide = write_molecule_slide(
        view,
        slide_dir / f"astro_molecules_{view_choice}.pptx",
        profile="balanced",
    )
    products.append(
        GeneratedProduct(
            label="ISM/CSM detections slide",
            path=ism_slide,
            description="Balanced PowerPoint slide of secure ISM/CSM molecules.",
        )
    )

    ism_layout = build_molecule_slide_layout(view, profile="balanced")
    ism_report = write_molecule_slide_report(
        ism_layout,
        slide_dir / f"astro_molecules_{view_choice}_layout.md",
    )
    products.append(
        GeneratedProduct(
            label="ISM/CSM slide layout report",
            path=ism_report,
            description="Text report for the generated ISM/CSM slide layout.",
        )
    )

    ppd_slide = write_ppd_detection_slide(
        view,
        slide_dir / f"ppd_molecules_{view_choice}.pptx",
    )
    products.append(
        GeneratedProduct(
            label="PPD detections slide",
            path=ppd_slide,
            description="PowerPoint slide of PPD molecules and isotopologues.",
        )
    )

    ppd_layout = build_ppd_detection_slide_layout(view)
    ppd_report = write_molecule_slide_report(
        ppd_layout,
        slide_dir / f"ppd_molecules_{view_choice}_layout.md",
    )
    products.append(
        GeneratedProduct(
            label="PPD slide layout report",
            path=ppd_report,
            description="Text report for the generated PPD slide layout.",
        )
    )
    return products


def _generate_tables(
    view: CensusView,
    *,
    output_dir: Path,
) -> list[GeneratedProduct]:
    """Generate standard LaTeX table fragments."""
    table_dir = output_dir / "tables"
    table_dir.mkdir(parents=True, exist_ok=True)
    writers = [
        ("Scalar LaTeX fragments", write_scalar_fragments),
        ("ISM/CSM molecule tables", lambda selected_view, out: write_ism_tables(selected_view, out, layout="balanced")),
        ("External-galaxy molecule table", write_exgal_table),
        ("PPD molecule table", write_ppd_table),
        ("Exoplanet-atmosphere molecule table", write_exoplanet_table),
        ("Interstellar-ice molecule table", write_ice_table),
        ("Detection-rate table", write_rate_by_atoms_table),
        ("Facility count table", write_facility_table),
        ("Source count table", write_source_table),
    ]
    products: list[GeneratedProduct] = []
    before = set(table_dir.glob("*.tex"))
    for label, writer in writers:
        writer(view, table_dir)
        after = set(table_dir.glob("*.tex"))
        for path in sorted(after - before):
            products.append(
                GeneratedProduct(
                    label=f"{label}: {path.name}",
                    path=path,
                    description="Generated LaTeX manuscript fragment.",
                )
            )
        before = after
    return products


def _write_index(
    output_dir: Path,
    products: list[GeneratedProduct],
    *,
    view_choice: str,
    db: Database,
) -> Path:
    """Write a simple static HTML index for GitHub Pages."""
    generated = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    groups = [
        ("Key Products", lambda path: path.suffix in {".zip", ".pptx"}),
        ("Figures", lambda path: "figures" in path.parts),
        ("Tables", lambda path: "tables" in path.parts),
        ("Reports", lambda path: path.suffix == ".md"),
    ]

    lines = [
        "<!doctype html>",
        '<html lang="en">',
        "<head>",
        '  <meta charset="utf-8">',
        '  <meta name="viewport" content="width=device-width, initial-scale=1">',
        "  <title>astromol latest outputs</title>",
        "  <style>",
        "    body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif; margin: 2rem auto; max-width: 980px; padding: 0 1rem; line-height: 1.45; }",
        "    h1, h2 { line-height: 1.15; }",
        "    .meta { color: #555; }",
        "    li { margin: 0.35rem 0; }",
        "    a { color: #005bbb; }",
        "  </style>",
        "</head>",
        "<body>",
        "  <h1>astromol latest generated outputs</h1>",
        f"  <p class=\"meta\">Generated {escape(generated)} from view <code>{escape(view_choice)}</code>.</p>",
        f"  <p class=\"meta\">Database loaded {len(db.molecules)} molecules, {len(db.detections)} detections, {len(db.sources)} sources, and {len(db.telescopes)} telescopes.</p>",
        "  <p>Use these files for the latest standard figures and slides. Use the Colab notebooks when you want custom views, alternate formats, or interactive regeneration.</p>",
    ]

    seen: set[Path] = set()
    for title, predicate in groups:
        selected = [
            product
            for product in products
            if product.path not in seen and predicate(product.path.relative_to(output_dir))
        ]
        if not selected:
            continue
        lines.append(f"  <h2>{escape(title)}</h2>")
        lines.append("  <ul>")
        for product in sorted(selected, key=lambda item: str(item.path)):
            seen.add(product.path)
            href = product.path.relative_to(output_dir).as_posix()
            lines.append(
                "    <li>"
                f"<a href=\"{escape(href)}\">{escape(product.label)}</a>"
                f" - {escape(product.description)}"
                "</li>"
            )
        lines.append("  </ul>")

    remaining = [product for product in products if product.path not in seen]
    if remaining:
        lines.append("  <h2>Other Outputs</h2>")
        lines.append("  <ul>")
        for product in sorted(remaining, key=lambda item: str(item.path)):
            href = product.path.relative_to(output_dir).as_posix()
            lines.append(
                "    <li>"
                f"<a href=\"{escape(href)}\">{escape(product.label)}</a>"
                f" - {escape(product.description)}"
                "</li>"
            )
        lines.append("  </ul>")

    lines.extend(["</body>", "</html>"])
    index_path = output_dir / "index.html"
    index_path.write_text("\n".join(lines) + "\n")
    return index_path


def _write_zip(output_dir: Path) -> Path:
    """Write a zip archive of generated outputs."""
    zip_path = output_dir / "astromol_latest_outputs.zip"
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(output_dir.rglob("*")):
            if path.is_file() and path != zip_path:
                archive.write(path, path.relative_to(output_dir))
    return zip_path


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments for the output generator."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory where generated products should be written.",
    )
    parser.add_argument(
        "--view",
        default="current",
        help="View to generate: current or a census year such as 2026.",
    )
    parser.add_argument(
        "--formats",
        nargs="+",
        default=list(DEFAULT_FORMATS),
        help="Figure formats to generate, usually png and/or pdf.",
    )
    parser.add_argument(
        "--no-clean",
        action="store_true",
        help="Do not remove the output directory before generating files.",
    )
    parser.add_argument(
        "--skip-tables",
        action="store_true",
        help="Skip LaTeX table-fragment generation.",
    )
    parser.add_argument(
        "--skip-slides",
        action="store_true",
        help="Skip PowerPoint slide generation.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """Command-line entry point for standard output generation."""
    args = parse_args(argv)
    products = generate_standard_outputs(
        args.output_dir,
        view_choice=args.view,
        formats=tuple(args.formats),
        clean=not args.no_clean,
        include_tables=not args.skip_tables,
        include_slides=not args.skip_slides,
    )
    print(f"Generated {len(products)} product(s) in {args.output_dir}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
