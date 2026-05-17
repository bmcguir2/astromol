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
