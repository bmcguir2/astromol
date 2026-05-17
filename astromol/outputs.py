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
    "other": "Other Files",
}

SECTION_DESCRIPTIONS = {
    "other": "Files that do not fit one of the main page sections.",
}


@dataclass(frozen=True)
class GeneratedProduct:
    """One generated output file for the static product index."""

    label: str
    path: Path
    description: str


@dataclass(frozen=True)
class FigureProductGroup:
    """One figure with its available download formats."""

    name: str
    label: str
    description: str
    png_product: GeneratedProduct | None = None
    pdf_product: GeneratedProduct | None = None


def _relative_product_path(output_dir: Path, product: GeneratedProduct) -> Path:
    """Return a product path relative to the bundle root."""
    return product.path.relative_to(output_dir)


def _product_section(relative_path: Path) -> str:
    """Return the inventory section key for one generated product."""
    if relative_path.name == "astromol_latest_outputs.zip":
        return "hidden"
    if relative_path.parts and relative_path.parts[0] == "figures":
        return "hidden"
    if relative_path.parts and relative_path.parts[0] == "slides":
        return "hidden"
    if relative_path.parts and relative_path.parts[0] == "tables":
        return "hidden"
    return "other"


def _featured_product_key(relative_path: Path) -> str | None:
    """Return the featured-product key for one path."""
    path_text = relative_path.as_posix()
    if relative_path.name == "astromol_latest_outputs.zip":
        return "bundle"
    if path_text.startswith("slides/astro_molecules_") and relative_path.suffix == ".pptx":
        return "ism_slide"
    if path_text.startswith("slides/ppd_molecules_") and relative_path.suffix == ".pptx":
        return "ppd_slide"
    if path_text == "figures/png/cumulative_detections.png":
        return "cumulative_figure"
    return None


def _featured_products(
    output_dir: Path,
    products: list[GeneratedProduct],
) -> dict[str, GeneratedProduct]:
    """Return the preferred top-of-page downloads when available."""
    selected: dict[str, GeneratedProduct] = {}
    for product in products:
        key = _featured_product_key(_relative_product_path(output_dir, product))
        if key is not None:
            selected[key] = product
    return selected


def _group_figure_products(
    output_dir: Path,
    products: list[GeneratedProduct],
) -> list[FigureProductGroup]:
    """Return figure downloads grouped by figure name across formats."""
    grouped: dict[str, dict[str, GeneratedProduct]] = {}
    for product in products:
        relative_path = _relative_product_path(output_dir, product)
        if not (relative_path.parts and relative_path.parts[0] == "figures"):
            continue
        group = grouped.setdefault(relative_path.stem, {})
        group[relative_path.suffix.lower()] = product

    figure_groups: list[FigureProductGroup] = []
    for name in sorted(grouped):
        formats = grouped[name]
        png_product = formats.get(".png")
        pdf_product = formats.get(".pdf")
        display_product = png_product or pdf_product
        if display_product is None:
            continue
        label = display_product.label.removesuffix(" (PNG)").removesuffix(" (PDF)")
        figure_groups.append(
            FigureProductGroup(
                name=name,
                label=label,
                description=display_product.description,
                png_product=png_product,
                pdf_product=pdf_product,
            )
        )
    return figure_groups


def _link_attributes(relative_path: Path) -> str:
    """Return extra anchor attributes for one generated product link."""
    if relative_path.suffix.lower() in {".png", ".jpg", ".jpeg", ".gif", ".svg", ".webp", ".pdf"}:
        return ' target="_blank" rel="noopener noreferrer"'
    return ""


def _render_link(output_dir: Path, product: GeneratedProduct, text: str, *, css_class: str = "") -> str:
    """Render one anchor element for a generated product."""
    relative_path = _relative_product_path(output_dir, product)
    class_attr = f' class="{css_class}"' if css_class else ""
    return (
        f'<a href="{escape(relative_path.as_posix())}"{class_attr}{_link_attributes(relative_path)}>'
        f"{escape(text)}</a>"
    )


def _render_product_item(output_dir: Path, product: GeneratedProduct) -> str:
    """Render one generated product as an HTML list item."""
    relative_path = _relative_product_path(output_dir, product)
    href = relative_path.as_posix()
    link_attrs = _link_attributes(relative_path)
    return (
        "        <li class=\"product-item\">"
        f"<a href=\"{escape(href)}\"{link_attrs}>{escape(product.label)}</a>"
        f"<span class=\"product-desc\">{escape(product.description)}</span>"
        f"<span class=\"product-path\">{escape(href)}</span>"
        "</li>"
    )


def _render_figure_card(output_dir: Path, group: FigureProductGroup) -> list[str]:
    """Render one figure card with PNG/PDF links."""
    lines = [
        "        <article class=\"figure-card\">",
        f"          <h3>{escape(group.label)}</h3>",
        f"          <p>{escape(group.description)}</p>",
        "          <div class=\"figure-links\">",
    ]
    if group.png_product is not None:
        lines.append(
            "            "
            + _render_link(output_dir, group.png_product, "PNG Preview", css_class="figure-link")
        )
    if group.pdf_product is not None:
        lines.append(
            "            "
            + _render_link(output_dir, group.pdf_product, "PDF File", css_class="figure-link")
        )
    lines.extend(
        [
            "          </div>",
            "        </article>",
        ]
    )
    return lines


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
    figure_groups = _group_figure_products(output_dir, products)
    section_order = ("other",)
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
        "figures": len(figure_groups),
        "slides": sum(
            1
            for product in products
            if _featured_product_key(_relative_product_path(output_dir, product)) in {"ism_slide", "ppd_slide"}
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
        "    .featured-card .tag, .bundle-card .tag, .slides-card .tag, .feature-card .tag { display: inline-block; margin-bottom: 10px; border-radius: 999px; padding: 5px 10px; background: rgba(24, 116, 208, 0.10); color: var(--blue-dark); font-size: 0.78rem; font-weight: 700; letter-spacing: 0.05em; text-transform: uppercase; }",
        "    .featured-card h3 { margin-bottom: 10px; }",
        "    .featured-card p { margin: 0 0 10px; color: var(--muted); }",
        "    .featured-card .file { display: block; color: var(--muted); font-size: 0.88rem; word-break: break-word; }",
        "    .primary-grid { display: grid; gap: 16px; grid-template-columns: repeat(3, minmax(0, 1fr)); align-items: start; }",
        "    .bundle-card, .slides-card, .feature-card, .figure-card, .inventory-card { background: var(--panel-strong); border: 1px solid var(--line); border-radius: 22px; box-shadow: var(--shadow); }",
        "    .bundle-card, .slides-card, .feature-card { padding: 18px 18px 16px; }",
        "    .slides-grid { display: grid; gap: 14px; grid-template-columns: 1fr; }",
        "    .slide-card { border: 1px solid rgba(17, 36, 58, 0.09); border-radius: 16px; padding: 14px; background: rgba(255, 255, 255, 0.55); }",
        "    .slide-card h3, .bundle-card h3, .feature-card h3 { margin-bottom: 8px; }",
        "    .slide-card p, .bundle-card p, .feature-card p { margin: 0 0 10px; color: var(--muted); }",
        "    .primary-link, .figure-link { display: inline-flex; align-items: center; justify-content: center; border-radius: 999px; padding: 9px 14px; border: 1px solid rgba(17, 36, 58, 0.12); background: rgba(24, 116, 208, 0.08); font-weight: 600; }",
        "    .primary-link:hover, .figure-link:hover { text-decoration: none; background: rgba(24, 116, 208, 0.14); }",
        "    .figure-grid { display: grid; gap: 16px; grid-template-columns: repeat(3, minmax(0, 1fr)); }",
        "    .figure-card { padding: 18px 18px 16px; }",
        "    .figure-card h3 { margin-bottom: 8px; }",
        "    .figure-card p { margin: 0 0 12px; color: var(--muted); }",
        "    .figure-links { display: flex; flex-wrap: wrap; gap: 10px; }",
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
        "    @media (max-width: 960px) { .figure-grid { grid-template-columns: repeat(2, minmax(0, 1fr)); } }",
        "    @media (max-width: 760px) { .page { padding: 22px 14px 38px; } .hero { padding: 22px 20px 20px; border-radius: 22px; } .hero-grid, .primary-grid, .slides-grid, .figure-grid { grid-template-columns: 1fr; } }",
        "  </style>",
        "</head>",
        "<body>",
        "  <main class=\"page\">",
        "    <header class=\"hero\">",
        "      <p class=\"eyebrow\">astromol generated outputs</p>",
        "      <h1>Latest standard figures, tables, and slides</h1>",
        "      <p class=\"lede\">The latest set of standard downloads generated from the selected census view. Start with the main files below, then browse the full inventory if you need another format or manuscript fragment.</p>",
        "      <div class=\"stat-row\">",
        f"        <span class=\"stat\">Generated <strong>{escape(generated)}</strong></span>",
        f"        <span class=\"stat\">View <code>{escape(view_choice)}</code></span>",
        f"        <span class=\"stat\">{product_counts['figures']} figures</span>",
        f"        <span class=\"stat\">{product_counts['slides']} slide deck(s)</span>",
        "      </div>",
        "      <div class=\"hero-grid\">",
        "        <section class=\"hero-panel\">",
        "          <h2>How To Use This Page</h2>",
        "          <p>Download the bundle if you want everything at once, open the slide decks for ready-to-use summaries, and browse the PNG or PDF sections when you only need a specific figure. Use the notebooks for custom views or selective regeneration.</p>",
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
    bundle_product = featured_products.get("bundle")
    ism_slide = featured_products.get("ism_slide")
    ppd_slide = featured_products.get("ppd_slide")
    cumulative_figure = next((group for group in figure_groups if group.name == "cumulative_detections"), None)
    public_figure_groups = [group for group in figure_groups if group.name != "cumulative_detections"]
    if bundle_product or cumulative_figure or ism_slide or ppd_slide:
        lines.extend(
            [
                "    <section class=\"section\">",
                "      <h2>Primary Downloads</h2>",
                "      <p class=\"section-copy\">Start with the full bundle, the standard cumulative figure, or the two standard slide decks.</p>",
                "      <div class=\"primary-grid\">",
            ]
        )
        if bundle_product is not None:
            lines.extend(
                [
                    "        <article class=\"bundle-card\">",
                    "          <span class=\"tag\">Featured</span>",
                    f"          <h3>{escape(bundle_product.label)}</h3>",
                    f"          <p>{escape(bundle_product.description)}</p>",
                    "          "
                    + _render_link(output_dir, bundle_product, "Download bundle", css_class="primary-link"),
                    "        </article>",
                ]
            )
        if cumulative_figure is not None:
            lines.extend(
                [
                    "        <article class=\"feature-card\">",
                    "          <span class=\"tag\">Figure</span>",
                    f"          <h3>{escape(cumulative_figure.label)}</h3>",
                    f"          <p>{escape(cumulative_figure.description)}</p>",
                    "          <div class=\"figure-links\">",
                ]
            )
            if cumulative_figure.png_product is not None:
                lines.append(
                    "            "
                    + _render_link(output_dir, cumulative_figure.png_product, "PNG Preview", css_class="figure-link")
                )
            if cumulative_figure.pdf_product is not None:
                lines.append(
                    "            "
                    + _render_link(output_dir, cumulative_figure.pdf_product, "PDF File", css_class="figure-link")
                )
            lines.extend(
                [
                    "          </div>",
                    "        </article>",
                ]
            )
        lines.extend(
            [
                "        <section class=\"slides-card\">",
                "          <span class=\"tag\">Slides</span>",
                "          <h3>Presentation Decks</h3>",
                "          <div class=\"slides-grid\">",
            ]
        )
        for product in (ism_slide, ppd_slide):
            if product is None:
                continue
            lines.extend(
                [
                    "            <article class=\"slide-card\">",
                    f"              <h3>{escape(product.label)}</h3>",
                    f"              <p>{escape(product.description)}</p>",
                    "              "
                    + _render_link(output_dir, product, "Open slide deck", css_class="primary-link"),
                    "            </article>",
                ]
            )
        lines.extend(
            [
                "          </div>",
                "        </section>",
                "      </div>",
                "    </section>",
            ]
        )

    if public_figure_groups:
        lines.extend(
            [
                "    <section class=\"section\">",
                "      <h2>Figures</h2>",
                "      <p class=\"section-copy\">Each figure card includes a PNG Preview link and, when available, a PDF File version of the same plot.</p>",
                "      <div class=\"figure-grid\">",
            ]
        )
        for group in public_figure_groups:
            lines.extend(_render_figure_card(output_dir, group))
        lines.extend(
            [
                "      </div>",
                "    </section>",
            ]
        )

    if any(section_products[section_key] for section_key in section_order):
        lines.extend(
            [
                "    <section class=\"section\">",
                "      <h2>Additional Files</h2>",
                "      <p class=\"section-copy\">Layout diagnostics and any uncategorized files are listed here.</p>",
                "      <div class=\"inventory-grid\">",
            ]
        )
        for section_key in section_order:
            lines.extend(_render_inventory_section(output_dir, section_key, section_products[section_key]))
        lines.extend(
            [
                "      </div>",
                "    </section>",
            ]
        )
    lines.extend(
        [
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
