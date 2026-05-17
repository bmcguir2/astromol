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
from .registry import (
    FIGURE_OUTPUTS,
    SLIDE_OUTPUTS,
    TABLE_OUTPUTS,
    OutputContext,
)


DEFAULT_OUTPUT_DIR = Path("build") / "astromol_outputs"
DEFAULT_FORMATS = ("png", "pdf")

SECTION_TITLES = {
    "bundle": "Bundle",
    "slides": "Slide Decks",
    "figures_png": "PNG Figures",
    "figures_pdf": "PDF Figures",
    "tables": "LaTeX Tables",
    "reports": "Layout Reports",
    "other": "Other Outputs",
}

SECTION_DESCRIPTIONS = {
    "bundle": "One archive containing the complete standard-output set.",
    "slides": "PowerPoint decks for the current standard presentation products.",
    "figures_png": "Quick-look raster figures for browsing and slide reuse.",
    "figures_pdf": "Publication-friendly vector figure files.",
    "tables": "Generated LaTeX fragments for manuscript tables and scalar inputs.",
    "reports": "Plain-text slide-layout diagnostics written alongside the decks.",
    "other": "Files that do not fit one of the main output groups.",
}


@dataclass(frozen=True)
class GeneratedProduct:
    """One generated output file for the static product index."""

    label: str
    path: Path
    description: str


def _relative_product_path(output_dir: Path, product: GeneratedProduct) -> Path:
    """Return a product path relative to the bundle root."""
    return product.path.relative_to(output_dir)


def _product_section(relative_path: Path) -> str:
    """Return the inventory section key for one generated product."""
    if relative_path.name == "astromol_latest_outputs.zip":
        return "bundle"
    if relative_path.parts and relative_path.parts[0] == "slides":
        return "reports" if relative_path.suffix == ".md" else "slides"
    if relative_path.parts and relative_path.parts[0] == "figures":
        return "figures_png" if relative_path.suffix == ".png" else "figures_pdf"
    if relative_path.parts and relative_path.parts[0] == "tables":
        return "tables"
    if relative_path.suffix == ".md":
        return "reports"
    return "other"


def _featured_product_key(relative_path: Path) -> tuple[str, int] | None:
    """Return the featured-product key and preference rank for one path."""
    path_text = relative_path.as_posix()
    if relative_path.name == "astromol_latest_outputs.zip":
        return ("bundle", 0)
    if path_text.startswith("slides/astro_molecules_") and relative_path.suffix == ".pptx":
        return ("ism_slide", 0)
    if path_text.startswith("slides/ppd_molecules_") and relative_path.suffix == ".pptx":
        return ("ppd_slide", 0)
    if path_text == "figures/png/cumulative_detections.png":
        return ("cumulative_figure", 0)
    if path_text == "figures/pdf/cumulative_detections.pdf":
        return ("cumulative_figure", 1)
    return None


def _featured_products(
    output_dir: Path,
    products: list[GeneratedProduct],
) -> list[GeneratedProduct]:
    """Return the preferred top-of-page downloads when available."""
    selected: dict[str, tuple[int, GeneratedProduct]] = {}
    for product in products:
        feature = _featured_product_key(_relative_product_path(output_dir, product))
        if feature is None:
            continue
        key, rank = feature
        current = selected.get(key)
        if current is None or rank < current[0]:
            selected[key] = (rank, product)

    order = ("bundle", "cumulative_figure", "ism_slide", "ppd_slide")
    return [selected[key][1] for key in order if key in selected]


def _render_product_item(output_dir: Path, product: GeneratedProduct) -> str:
    """Render one generated product as an HTML list item."""
    href = _relative_product_path(output_dir, product).as_posix()
    return (
        "        <li class=\"product-item\">"
        f"<a href=\"{escape(href)}\">{escape(product.label)}</a>"
        f"<span class=\"product-desc\">{escape(product.description)}</span>"
        f"<span class=\"product-path\">{escape(href)}</span>"
        "</li>"
    )


def _render_inventory_section(
    output_dir: Path,
    section_key: str,
    products: list[GeneratedProduct],
) -> list[str]:
    """Render one inventory section card."""
    if not products:
        return []
    lines = [
        "      <section class=\"inventory-card\">",
        f"        <h3>{escape(SECTION_TITLES[section_key])}</h3>",
        f"        <p>{escape(SECTION_DESCRIPTIONS[section_key])}</p>",
        "        <ul class=\"product-list\">",
    ]
    lines.extend(_render_product_item(output_dir, product) for product in products)
    lines.extend(
        [
            "        </ul>",
            "      </section>",
        ]
    )
    return lines


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
    context = OutputContext(
        view=view,
        view_choice=view_choice,
        baseline_view=view_2021,
    )
    products: list[GeneratedProduct] = []

    products.extend(
        _generate_figures(
            context,
            output_dir=output_dir,
            formats=formats,
        )
    )
    if include_slides:
        products.extend(
            _generate_slides(
                context,
                output_dir=output_dir,
            )
        )
    if include_tables:
        products.extend(
            _generate_tables(
                context,
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
    context: OutputContext,
    *,
    output_dir: Path,
    formats: tuple[str, ...],
) -> list[GeneratedProduct]:
    """Generate the standard manuscript-facing figure set."""
    figure_dir = output_dir / "figures"
    products: list[GeneratedProduct] = []

    for spec in FIGURE_OUTPUTS:
        for file_format in formats:
            path = figure_dir / file_format / f"{spec.name}.{file_format}"
            path.parent.mkdir(parents=True, exist_ok=True)
            spec.write(context, path)
            products.append(
                GeneratedProduct(
                    label=f"{spec.label} ({file_format.upper()})",
                    path=path,
                    description=spec.description,
                )
            )
    return products


def _generate_slides(
    context: OutputContext,
    *,
    output_dir: Path,
) -> list[GeneratedProduct]:
    """Generate standard PowerPoint slides and layout reports."""
    slide_dir = output_dir / "slides"
    slide_dir.mkdir(parents=True, exist_ok=True)
    products: list[GeneratedProduct] = []

    for spec in SLIDE_OUTPUTS:
        slide_path = spec.write_slide(context, slide_dir)
        products.append(
            GeneratedProduct(
                label=spec.label,
                path=slide_path,
                description=spec.description,
            )
        )

        report_path = spec.write_report(context, slide_dir)
        if report_path is not None:
            products.append(
                GeneratedProduct(
                    label=spec.report_label or f"{spec.label} layout report",
                    path=report_path,
                    description=spec.report_description or "Text report for the generated slide layout.",
                )
            )
    return products


def _generate_tables(
    context: OutputContext,
    *,
    output_dir: Path,
) -> list[GeneratedProduct]:
    """Generate standard LaTeX table fragments."""
    table_dir = output_dir / "tables"
    table_dir.mkdir(parents=True, exist_ok=True)
    products: list[GeneratedProduct] = []
    for spec in TABLE_OUTPUTS:
        for path in spec.write(context.view, table_dir):
            products.append(
                GeneratedProduct(
                    label=f"{spec.label}: {path.name}",
                    path=path,
                    description=spec.description,
                )
            )
    return products


def _write_index(
    output_dir: Path,
    products: list[GeneratedProduct],
    *,
    view_choice: str,
    db: Database,
) -> Path:
    """Write the static GitHub Pages landing page for generated outputs."""
    generated = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    featured_products = _featured_products(output_dir, products)
    section_order = (
        "bundle",
        "slides",
        "figures_png",
        "figures_pdf",
        "tables",
        "reports",
        "other",
    )
    section_products = {
        key: sorted(
            [
                product
                for product in products
                if _product_section(_relative_product_path(output_dir, product)) == key
            ],
            key=lambda item: str(_relative_product_path(output_dir, item)),
        )
        for key in section_order
    }
    product_counts = {
        "figures": sum(
            1
            for product in products
            if _product_section(_relative_product_path(output_dir, product))
            in {"figures_png", "figures_pdf"}
        ),
        "slides": sum(
            1
            for product in products
            if _product_section(_relative_product_path(output_dir, product)) == "slides"
        ),
        "tables": sum(
            1
            for product in products
            if _product_section(_relative_product_path(output_dir, product)) == "tables"
        ),
        "reports": sum(
            1
            for product in products
            if _product_section(_relative_product_path(output_dir, product)) == "reports"
        ),
    }

    lines = [
        "<!doctype html>",
        '<html lang="en">',
        "<head>",
        '  <meta charset="utf-8">',
        '  <meta name="viewport" content="width=device-width, initial-scale=1">',
        "  <title>astromol latest outputs</title>",
        "  <style>",
        "    :root { color-scheme: light; --bg: #f4f0e8; --panel: rgba(255, 253, 248, 0.92); --panel-strong: rgba(255, 252, 246, 0.98); --ink: #11243a; --muted: #53606d; --line: rgba(17, 36, 58, 0.14); --blue: #1874d0; --blue-dark: #0e4a85; --accent: #8f2d2d; --shadow: 0 18px 46px rgba(17, 36, 58, 0.10); }",
        "    * { box-sizing: border-box; }",
        "    body { margin: 0; background: radial-gradient(circle at top left, rgba(24, 116, 208, 0.14), transparent 34%), linear-gradient(180deg, #f8f4ed 0%, #f1ece3 48%, #ebe6dc 100%); color: var(--ink); font-family: 'Avenir Next', 'Segoe UI', sans-serif; line-height: 1.55; }",
        "    a { color: var(--blue-dark); text-decoration: none; }",
        "    a:hover { text-decoration: underline; }",
        "    code { font-family: 'SFMono-Regular', Consolas, monospace; font-size: 0.96em; }",
        "    .page { max-width: 1120px; margin: 0 auto; padding: 32px 18px 48px; }",
        "    .hero { background: linear-gradient(140deg, rgba(255, 252, 246, 0.97), rgba(246, 248, 252, 0.92)); border: 1px solid var(--line); border-radius: 28px; box-shadow: var(--shadow); padding: 28px 30px 24px; }",
        "    .eyebrow { margin: 0 0 10px; color: var(--accent); font-size: 0.82rem; font-weight: 700; letter-spacing: 0.12em; text-transform: uppercase; }",
        "    h1, h2, h3 { font-family: Georgia, 'Iowan Old Style', serif; line-height: 1.08; letter-spacing: -0.01em; }",
        "    h1 { margin: 0; font-size: clamp(2.2rem, 4vw, 3.35rem); max-width: 13ch; }",
        "    h2 { margin: 0 0 10px; font-size: 1.55rem; }",
        "    h3 { margin: 0 0 8px; font-size: 1.18rem; }",
        "    .lede { max-width: 70ch; margin: 14px 0 0; color: var(--muted); font-size: 1.03rem; }",
        "    .hero-grid { display: grid; gap: 18px; grid-template-columns: minmax(0, 1.6fr) minmax(280px, 1fr); margin-top: 22px; }",
        "    .hero-panel { background: var(--panel); border: 1px solid var(--line); border-radius: 20px; padding: 18px 20px; }",
        "    .hero-panel h2 { font-size: 1.1rem; margin-bottom: 8px; }",
        "    .hero-panel p { margin: 0; color: var(--muted); }",
        "    .fact-list { display: grid; gap: 10px; margin: 0; }",
        "    .fact { display: flex; justify-content: space-between; gap: 12px; border-bottom: 1px solid rgba(17, 36, 58, 0.08); padding-bottom: 10px; }",
        "    .fact:last-child { border-bottom: none; padding-bottom: 0; }",
        "    .fact span:first-child { color: var(--muted); }",
        "    .fact strong { font-weight: 600; text-align: right; }",
        "    .stat-row { display: flex; flex-wrap: wrap; gap: 10px; margin-top: 18px; }",
        "    .stat { border: 1px solid var(--line); border-radius: 999px; background: rgba(255, 255, 255, 0.74); padding: 8px 12px; font-size: 0.94rem; }",
        "    .section { margin-top: 28px; }",
        "    .section-copy { max-width: 74ch; margin: 0 0 16px; color: var(--muted); }",
        "    .featured-grid { display: grid; gap: 16px; grid-template-columns: repeat(auto-fit, minmax(220px, 1fr)); }",
        "    .featured-card, .inventory-card { background: var(--panel-strong); border: 1px solid var(--line); border-radius: 22px; box-shadow: var(--shadow); }",
        "    .featured-card { padding: 18px 18px 16px; }",
        "    .featured-card .tag { display: inline-block; margin-bottom: 10px; border-radius: 999px; padding: 5px 10px; background: rgba(24, 116, 208, 0.10); color: var(--blue-dark); font-size: 0.78rem; font-weight: 700; letter-spacing: 0.05em; text-transform: uppercase; }",
        "    .featured-card h3 { margin-bottom: 10px; }",
        "    .featured-card p { margin: 0 0 10px; color: var(--muted); }",
        "    .featured-card .file { display: block; color: var(--muted); font-size: 0.88rem; word-break: break-word; }",
        "    .inventory-grid { display: grid; gap: 16px; grid-template-columns: repeat(auto-fit, minmax(280px, 1fr)); align-items: start; }",
        "    .inventory-card { padding: 18px 18px 16px; }",
        "    .inventory-card p { margin: 0 0 12px; color: var(--muted); }",
        "    .product-list { list-style: none; margin: 0; padding: 0; }",
        "    .product-item { padding: 10px 0; border-top: 1px solid rgba(17, 36, 58, 0.09); }",
        "    .product-item:first-child { border-top: none; padding-top: 0; }",
        "    .product-item a { font-weight: 600; }",
        "    .product-desc, .product-path { display: block; }",
        "    .product-desc { margin-top: 4px; color: var(--muted); }",
        "    .product-path { margin-top: 5px; font-size: 0.85rem; color: var(--muted); word-break: break-word; }",
        "    @media (max-width: 760px) { .page { padding: 22px 14px 38px; } .hero { padding: 22px 20px 20px; border-radius: 22px; } .hero-grid { grid-template-columns: 1fr; } }",
        "  </style>",
        "</head>",
        "<body>",
        "  <main class=\"page\">",
        "    <header class=\"hero\">",
        "      <p class=\"eyebrow\">astromol generated outputs</p>",
        "      <h1>Latest standard figures, tables, and slides</h1>",
        "      <p class=\"lede\">Registry-driven standard bundle for the selected census view. Start with the primary downloads below, then use the full inventory for alternate formats and manuscript fragments.</p>",
        "      <div class=\"stat-row\">",
        f"        <span class=\"stat\">Generated <strong>{escape(generated)}</strong></span>",
        f"        <span class=\"stat\">View <code>{escape(view_choice)}</code></span>",
        f"        <span class=\"stat\">{product_counts['figures']} figure files</span>",
        f"        <span class=\"stat\">{product_counts['slides']} slide deck(s)</span>",
        f"        <span class=\"stat\">{product_counts['tables']} table fragment(s)</span>",
        f"        <span class=\"stat\">{product_counts['reports']} report(s)</span>",
        "      </div>",
        "      <div class=\"hero-grid\">",
        "        <section class=\"hero-panel\">",
        "          <h2>How To Use This Page</h2>",
        "          <p>Download the bundle for everything at once, open the slide decks for presentation-ready summaries, and browse the PNG/PDF figure sections when you need one product in a specific format. Use the notebooks for custom views or selective regeneration.</p>",
        "        </section>",
        "        <section class=\"hero-panel\">",
        "          <h2>Data Snapshot</h2>",
        "          <div class=\"fact-list\">",
        f"            <div class=\"fact\"><span>Molecules</span><strong>{len(db.molecules)}</strong></div>",
        f"            <div class=\"fact\"><span>Detections</span><strong>{len(db.detections)}</strong></div>",
        f"            <div class=\"fact\"><span>Sources</span><strong>{len(db.sources)}</strong></div>",
        f"            <div class=\"fact\"><span>Telescopes</span><strong>{len(db.telescopes)}</strong></div>",
        "          </div>",
        "        </section>",
        "      </div>",
        "    </header>",
    ]
    if featured_products:
        lines.extend(
            [
                "    <section class=\"section\">",
                "      <h2>Primary Downloads</h2>",
                "      <p class=\"section-copy\">The standard bundle and the most commonly requested presentation products are linked here first for quick access.</p>",
                "      <div class=\"featured-grid\">",
            ]
        )
        for product in featured_products:
            href = _relative_product_path(output_dir, product).as_posix()
            lines.extend(
                [
                    "        <article class=\"featured-card\">",
                    "          <span class=\"tag\">Featured</span>",
                    f"          <h3><a href=\"{escape(href)}\">{escape(product.label)}</a></h3>",
                    f"          <p>{escape(product.description)}</p>",
                    f"          <span class=\"file\">{escape(href)}</span>",
                    "        </article>",
                ]
            )
        lines.extend(
            [
                "      </div>",
                "    </section>",
            ]
        )

    lines.extend(
        [
            "    <section class=\"section\">",
            "      <h2>Full Inventory</h2>",
            "      <p class=\"section-copy\">Every file in the registry-driven bundle is grouped below by product type so alternate formats, manuscript fragments, and slide diagnostics are easy to find.</p>",
            "      <div class=\"inventory-grid\">",
        ]
    )
    for section_key in section_order:
        lines.extend(_render_inventory_section(output_dir, section_key, section_products[section_key]))
    lines.extend(
        [
            "      </div>",
            "    </section>",
            "  </main>",
            "</body>",
            "</html>",
        ]
    )
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
