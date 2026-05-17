"""Public command-line interface for astromol outputs."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLBACKEND", "Agg")

from .census import CensusView
from .database import Database
from .outputs import DEFAULT_FORMATS, DEFAULT_OUTPUT_DIR, generate_standard_outputs
from .registry import FIGURE_OUTPUTS, SLIDE_OUTPUTS, TABLE_OUTPUTS, OutputContext


def _names(specs: tuple[Any, ...]) -> tuple[str, ...]:
    """Return registry names for an argparse choices list."""
    return tuple(spec.name for spec in specs)


def _spec_by_name(specs: tuple[Any, ...], name: str) -> Any:
    """Return one registry spec by stable name."""
    for spec in specs:
        if spec.name == name:
            return spec
    valid_names = ", ".join(_names(specs))
    raise ValueError(f"Unknown output name {name!r}. Valid names: {valid_names}")


def _build_context(view_choice: str) -> OutputContext:
    """Build an output context for one CLI invocation."""
    db = Database()
    if view_choice == "current":
        view = CensusView.current(db)
    else:
        view = CensusView.for_census(db, view_choice)
    return OutputContext(
        view=view,
        view_choice=view_choice,
        baseline_view=CensusView.for_census(db, "2021"),
    )


def write_figure(name: str, output_path: str | Path, *, view_choice: str = "current") -> Path:
    """Write one standard figure selected by registry name."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    context = _build_context(view_choice)
    spec = _spec_by_name(FIGURE_OUTPUTS, name)
    return spec.write(context, output_path)


def write_table(name: str, output_dir: str | Path, *, view_choice: str = "current") -> list[Path]:
    """Write one standard table-fragment group selected by registry name."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    context = _build_context(view_choice)
    spec = _spec_by_name(TABLE_OUTPUTS, name)
    return spec.write(context.view, output_dir)


def write_slide(
    name: str,
    output_dir: str | Path,
    *,
    view_choice: str = "current",
    include_report: bool = False,
) -> list[Path]:
    """Write one standard slide deck selected by registry name."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    context = _build_context(view_choice)
    spec = _spec_by_name(SLIDE_OUTPUTS, name)

    paths = [spec.write_slide(context, output_dir)]
    if include_report:
        report_path = spec.write_report(context, output_dir)
        if report_path is not None:
            paths.append(report_path)
    return paths


def _print_registry(category: str) -> None:
    """Print standard output names grouped by registry category."""
    groups = []
    if category in {"all", "figures"}:
        groups.append(("Figures", FIGURE_OUTPUTS))
    if category in {"all", "tables"}:
        groups.append(("Tables", TABLE_OUTPUTS))
    if category in {"all", "slides"}:
        groups.append(("Slides", SLIDE_OUTPUTS))

    for index, (title, specs) in enumerate(groups):
        if index:
            print()
        print(f"{title}:")
        for spec in specs:
            print(f"  {spec.name}")
            print(f"    {spec.label}")
            print(f"    {spec.description}")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse public astromol CLI arguments."""
    parser = argparse.ArgumentParser(prog="astromol")
    subparsers = parser.add_subparsers(dest="command", required=True)

    list_parser = subparsers.add_parser("list", help="List standard output names.")
    list_parser.add_argument(
        "category",
        nargs="?",
        default="all",
        choices=("all", "figures", "tables", "slides"),
        help="Output category to list.",
    )

    figure_parser = subparsers.add_parser("figure", help="Generate one standard figure.")
    figure_parser.add_argument("name", choices=_names(FIGURE_OUTPUTS), help="Standard figure name.")
    figure_parser.add_argument(
        "--output",
        "-o",
        type=Path,
        required=True,
        help="Output figure path. The suffix selects the file format, such as .pdf or .png.",
    )
    figure_parser.add_argument(
        "--view",
        default="current",
        help="View to generate: current or a census year such as 2026.",
    )

    table_parser = subparsers.add_parser("table", help="Generate one standard table group.")
    table_parser.add_argument("name", choices=_names(TABLE_OUTPUTS), help="Standard table name.")
    table_parser.add_argument(
        "--output-dir",
        "-o",
        type=Path,
        default=Path("."),
        help="Directory where generated .tex fragments should be written.",
    )
    table_parser.add_argument(
        "--view",
        default="current",
        help="View to generate: current or a census year such as 2026.",
    )

    slide_parser = subparsers.add_parser("slide", help="Generate one standard slide deck.")
    slide_parser.add_argument("name", choices=_names(SLIDE_OUTPUTS), help="Standard slide name.")
    slide_parser.add_argument(
        "--output-dir",
        "-o",
        type=Path,
        default=Path("."),
        help="Directory where the generated PowerPoint file should be written.",
    )
    slide_parser.add_argument(
        "--view",
        default="current",
        help="View to generate: current or a census year such as 2026.",
    )
    slide_parser.add_argument(
        "--report",
        action="store_true",
        help="Also write the slide layout report when available.",
    )

    outputs_parser = subparsers.add_parser(
        "outputs",
        help="Generate the complete standard output bundle.",
    )
    outputs_parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory where generated products should be written.",
    )
    outputs_parser.add_argument(
        "--view",
        default="current",
        help="View to generate: current or a census year such as 2026.",
    )
    outputs_parser.add_argument(
        "--formats",
        nargs="+",
        default=list(DEFAULT_FORMATS),
        help="Figure formats to generate, usually png and/or pdf.",
    )
    outputs_parser.add_argument(
        "--no-clean",
        action="store_true",
        help="Do not remove the output directory before generating files.",
    )
    outputs_parser.add_argument(
        "--skip-tables",
        action="store_true",
        help="Skip LaTeX table-fragment generation.",
    )
    outputs_parser.add_argument(
        "--skip-slides",
        action="store_true",
        help="Skip PowerPoint slide generation.",
    )

    args = parser.parse_args(argv)
    if args.command == "figure" and not args.output.suffix:
        parser.error("figure --output must include a file suffix such as .pdf or .png")
    return args


def main(argv: list[str] | None = None) -> int:
    """Run the public astromol CLI."""
    args = parse_args(argv)

    if args.command == "list":
        _print_registry(args.category)
        return 0

    if args.command == "figure":
        path = write_figure(args.name, args.output, view_choice=args.view)
        print(path)
        return 0

    if args.command == "table":
        paths = write_table(args.name, args.output_dir, view_choice=args.view)
        for path in paths:
            print(path)
        return 0

    if args.command == "slide":
        paths = write_slide(
            args.name,
            args.output_dir,
            view_choice=args.view,
            include_report=args.report,
        )
        for path in paths:
            print(path)
        return 0

    if args.command == "outputs":
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

    raise ValueError(f"Unhandled command: {args.command}")


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
