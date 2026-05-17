"""Registry of standard astromol generated outputs.

The registry is intentionally metadata-first: each entry has a stable name,
human-readable label, description, and the callable pieces needed to generate
the product. This gives docs, notebooks, automated bundles, and future CLI
commands one canonical list of standard outputs.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .census import CensusView
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


DataBuilder = Callable[["OutputContext"], Any]
KeywordBuilder = Callable[["OutputContext"], Mapping[str, Any]]
FigureWriter = Callable[..., Path]
TableWriter = Callable[[CensusView, Path], Any]
SlideWriter = Callable[..., Path]
LayoutBuilder = Callable[..., Any]


def _empty_kwargs(_: "OutputContext") -> Mapping[str, Any]:
    """Return no keyword arguments for a registry entry."""
    return {}


@dataclass
class OutputContext:
    """Shared context for generating standard outputs from one selected view."""

    view: CensusView
    view_choice: str = "current"
    baseline_view: CensusView | None = None
    cache: dict[str, Any] = field(default_factory=dict)

    def cached(self, key: str, builder: Callable[[], Any]) -> Any:
        """Return a cached intermediate data product."""
        if key not in self.cache:
            self.cache[key] = builder()
        return self.cache[key]


@dataclass(frozen=True)
class FigureOutputSpec:
    """One standard figure product."""

    name: str
    label: str
    description: str
    data_builder: DataBuilder
    writer: FigureWriter
    writer_kwargs: KeywordBuilder = _empty_kwargs
    baseline_data_builder: DataBuilder | None = None

    def build_data(self, context: OutputContext) -> Any:
        """Build the primary data object for this figure."""
        return self.data_builder(context)

    def write(self, context: OutputContext, output_path: str | Path) -> Path:
        """Write this figure to ``output_path``."""
        data = self.build_data(context)
        kwargs = dict(self.writer_kwargs(context))
        if self.baseline_data_builder is None:
            return self.writer(data, output_path, **kwargs)
        baseline_data = self.baseline_data_builder(context)
        return self.writer(data, baseline_data, output_path, **kwargs)


@dataclass(frozen=True)
class TableOutputSpec:
    """One standard LaTeX table-fragment product group."""

    name: str
    label: str
    description: str
    writer: TableWriter

    def write(self, view: CensusView, output_dir: str | Path) -> list[Path]:
        """Write this table group and return generated ``.tex`` paths."""
        output_dir = Path(output_dir)
        before = set(output_dir.glob("*.tex"))
        result = self.writer(view, output_dir)
        if isinstance(result, Mapping):
            return sorted(output_dir / filename for filename in result)
        if result is not None:
            return sorted(Path(path) for path in result)
        after = set(output_dir.glob("*.tex"))
        return sorted(after - before)


@dataclass(frozen=True)
class SlideOutputSpec:
    """One standard PowerPoint slide plus its optional layout report."""

    name: str
    label: str
    description: str
    filename_template: str
    writer: SlideWriter
    writer_kwargs: KeywordBuilder = _empty_kwargs
    layout_builder: LayoutBuilder | None = None
    report_label: str | None = None
    report_description: str | None = None
    report_filename_template: str | None = None

    def filename(self, context: OutputContext) -> str:
        """Return this slide's filename for the selected context."""
        return self.filename_template.format(view_choice=context.view_choice)

    def report_filename(self, context: OutputContext) -> str | None:
        """Return this slide's layout-report filename when configured."""
        if self.report_filename_template is None:
            return None
        return self.report_filename_template.format(view_choice=context.view_choice)

    def write_slide(self, context: OutputContext, output_dir: str | Path) -> Path:
        """Write this PowerPoint slide."""
        output_path = Path(output_dir) / self.filename(context)
        return self.writer(context.view, output_path, **dict(self.writer_kwargs(context)))

    def write_report(self, context: OutputContext, output_dir: str | Path) -> Path | None:
        """Write this slide's layout report when configured."""
        report_filename = self.report_filename(context)
        if self.layout_builder is None or report_filename is None:
            return None
        layout = self.layout_builder(context.view, **dict(self.writer_kwargs(context)))
        return write_molecule_slide_report(layout, Path(output_dir) / report_filename)


def _view_2021(context: OutputContext) -> CensusView:
    """Return the 2021 comparison view for registry entries that need it."""
    if context.baseline_view is None:
        raise ValueError("A 2021 baseline view is required for this output.")
    return context.baseline_view


def _by_atoms(context: OutputContext) -> Any:
    return context.cached("cumulative_by_atoms", lambda: cumulative_by_atoms_data(context.view))


def _du_source(context: OutputContext) -> Any:
    return context.cached("du_by_source_type", lambda: du_by_source_type_data(context.view))


def _relative_du_source(context: OutputContext) -> Any:
    return context.cached(
        "relative_du_by_source_type",
        lambda: relative_du_by_source_type_data(context.view),
    )


def _mass_source(context: OutputContext) -> Any:
    return context.cached("mass_by_source_type", lambda: mass_by_source_type_data(context.view))


def _wavelength_source(context: OutputContext) -> Any:
    return context.cached(
        "wavelength_by_source_type",
        lambda: wavelength_by_source_type_data(context.view),
    )


def _rate_data(context: OutputContext) -> Any:
    return context.cached("detection_rate_by_atoms", lambda: detection_rate_by_atoms_data(context.view))


def _rate_2021(context: OutputContext) -> Any:
    return context.cached(
        "detection_rate_by_atoms_2021",
        lambda: detection_rate_by_atoms_data(_view_2021(context)),
    )


def _comparison_kwargs(context: OutputContext) -> Mapping[str, Any]:
    return {"current_label": context.view_choice}


def _modern_scopes_kwargs(_: OutputContext) -> Mapping[str, Any]:
    return {"style": "modern"}


def _balanced_slide_kwargs(_: OutputContext) -> Mapping[str, Any]:
    return {"profile": "balanced"}


def _write_ism_tables_balanced(view: CensusView, output_dir: Path) -> Any:
    return write_ism_tables(view, output_dir, layout="balanced")


FIGURE_OUTPUTS: tuple[FigureOutputSpec, ...] = (
    FigureOutputSpec(
        "cumulative_detections",
        "Cumulative ISM/CSM detections",
        "Secure, non-isotopologue ISM/CSM first detections through the selected view boundary.",
        lambda context: cumulative_detection_data(context.view),
        write_cumulative_detections_plot,
    ),
    FigureOutputSpec(
        "cumulative_by_atoms",
        "Cumulative detections by atom-count category",
        "Legacy-style line traces by atom-count category.",
        _by_atoms,
        write_cumulative_by_atoms_plot,
    ),
    FigureOutputSpec(
        "cumulative_by_atoms_stacked",
        "Stacked cumulative detections by atom-count category",
        "Production companion view showing how atom-count categories contribute to the cumulative inventory.",
        _by_atoms,
        write_stacked_cumulative_by_atoms_plot,
    ),
    FigureOutputSpec(
        "rolling_rate_by_atoms",
        "Rolling detection-rate heatmap by atom-count category",
        "Trailing rolling detections per year by atom-count category.",
        lambda context: rolling_rate_by_atoms_heatmap_data(context.view),
        write_rolling_rate_by_atoms_heatmap,
    ),
    FigureOutputSpec(
        "periodic_heatmap",
        "Periodic-table heatmap of elements in detected molecules",
        "Counts how many selected molecules contain each element.",
        lambda context: periodic_heatmap_data(context.view),
        write_periodic_heatmap,
    ),
    FigureOutputSpec(
        "du_bar_chart",
        "Degree-of-unsaturation exact-value bar chart",
        "Production exact-value degree-of-unsaturation bar chart.",
        lambda context: du_histogram_data(context.view),
        write_du_bar_chart,
    ),
    FigureOutputSpec(
        "kappa_distribution",
        "Ray asymmetry-parameter distribution",
        "Kappa distribution for molecules with usable rotational constants; linear molecules are included at kappa = -1.",
        lambda context: kappa_histogram_data(context.view),
        write_kappa_histogram,
    ),
    FigureOutputSpec(
        "molecule_type_pie",
        "Molecule-type pie chart",
        "Legacy categorical molecule-type summary.",
        lambda context: molecule_type_data(context.view),
        write_type_pie_chart,
    ),
    FigureOutputSpec(
        "source_type_pie",
        "Source-type pie chart",
        "Legacy generalized first-detection source-type summary.",
        lambda context: source_type_data(context.view),
        write_source_pie_chart,
    ),
    FigureOutputSpec(
        "individual_source_pie",
        "Individual-source contribution pie chart",
        "Legacy summary for major individual first-detection sources.",
        lambda context: individual_source_data(context.view),
        write_individual_source_pie_chart,
    ),
    FigureOutputSpec(
        "molecule_type_source_enrichment",
        "Molecule-type enrichment by source type",
        "Production enrichment matrix for molecule-type categories by generalized first-detection source type.",
        lambda context: molecule_type_by_source_type_data(context.view),
        write_molecule_type_by_source_enrichment_matrix,
    ),
    FigureOutputSpec(
        "du_by_source_boxplot",
        "Degree of unsaturation by source type",
        "Production source-category comparison for discrete degree-of-unsaturation values.",
        _du_source,
        write_du_by_source_type_boxplot,
    ),
    FigureOutputSpec(
        "relative_du_by_source_boxplot",
        "Relative degree of unsaturation by source type",
        "Production source-category comparison for relative degree-of-unsaturation values.",
        _relative_du_source,
        write_relative_du_by_source_type_boxplot,
    ),
    FigureOutputSpec(
        "mass_by_source_boxplot",
        "Molecular mass by source type",
        "Production source-category comparison for molecular masses.",
        _mass_source,
        write_mass_by_source_type_boxplot,
    ),
    FigureOutputSpec(
        "wavelength_by_source_stacked_bar",
        "First-detection wavelength share by source type",
        "Production replacement for the wavelength/source pie-grid.",
        _wavelength_source,
        write_wavelength_by_source_type_stacked_bar,
    ),
    FigureOutputSpec(
        "mass_by_wavelength_kde",
        "Molecular mass distribution by first-detection wavelength",
        "Legacy KDE comparison by first-detection wavelength.",
        lambda context: mass_by_wavelength_data(context.view, include_fullerenes=False),
        write_mass_by_wavelength_plot,
    ),
    FigureOutputSpec(
        "mass_by_wavelength_boxplot",
        "Molecular mass by first-detection wavelength",
        "Production wavelength comparison for molecular masses.",
        lambda context: mass_by_wavelength_data(context.view, include_fullerenes=False),
        write_mass_by_wavelength_boxplot,
    ),
    FigureOutputSpec(
        "wavelength_atoms_bubble_heatmap",
        "Detection wavelength by atom-count bubble heatmap",
        "Production replacement for the wavelength/atom-count KDE figure.",
        lambda context: molecules_by_wavelength_atoms_data(context.view),
        write_molecules_by_wavelength_atoms_bubble_heatmap,
    ),
    FigureOutputSpec(
        "detection_rate_by_atoms",
        "Detection rate by atom-count category",
        "Average detections per year by atom-count category.",
        _rate_data,
        write_detection_rate_by_atoms_plot,
    ),
    FigureOutputSpec(
        "detection_rate_by_atoms_comparison",
        "Detection rate by atom-count category with 2021 baseline",
        "Production comparison plot using the selected view as foreground and the 2021 census as background reference.",
        _rate_data,
        write_detection_rate_by_atoms_comparison_plot,
        writer_kwargs=_comparison_kwargs,
        baseline_data_builder=_rate_2021,
    ),
    FigureOutputSpec(
        "facility_share_bars",
        "Facility contribution share by facility era",
        "Production replacement for the facility-share pie chart.",
        lambda context: facility_share_data(context.view),
        write_facility_share_bars_plot,
    ),
    FigureOutputSpec(
        "scopes_by_year",
        "Cumulative facility contributions by year",
        "Modernized facility-contribution trace plot.",
        lambda context: scopes_by_year_data(context.view),
        write_scopes_by_year_plot,
        writer_kwargs=_modern_scopes_kwargs,
    ),
)


TABLE_OUTPUTS: tuple[TableOutputSpec, ...] = (
    TableOutputSpec(
        "scalars",
        "Scalar LaTeX fragments",
        "Individual scalar LaTeX inputs such as molecule counts, source counts, and detection rates.",
        write_scalar_fragments,
    ),
    TableOutputSpec(
        "ism_tables",
        "ISM/CSM molecule tables",
        "Balanced ISM/CSM molecule tables.",
        _write_ism_tables_balanced,
    ),
    TableOutputSpec(
        "exgal_table",
        "External-galaxy molecule table",
        "External-galaxy molecule table, isotopologue-free.",
        write_exgal_table,
    ),
    TableOutputSpec(
        "ppd_table",
        "PPD molecule table",
        "PPD molecule table including isotopologues by default.",
        write_ppd_table,
    ),
    TableOutputSpec(
        "exoplanet_table",
        "Exoplanet-atmosphere molecule table",
        "Exoplanet-atmosphere molecule table.",
        write_exoplanet_table,
    ),
    TableOutputSpec(
        "ice_table",
        "Interstellar-ice molecule table",
        "Interstellar-ice molecule table.",
        write_ice_table,
    ),
    TableOutputSpec(
        "rate_by_atoms_table",
        "Detection-rate table",
        "Detection-rate fit table by atom-count category.",
        write_rate_by_atoms_table,
    ),
    TableOutputSpec(
        "facility_table",
        "Facility count table",
        "Facility contribution count table.",
        write_facility_table,
    ),
    TableOutputSpec(
        "source_table",
        "Source count table",
        "Source contribution count table.",
        write_source_table,
    ),
)


SLIDE_OUTPUTS: tuple[SlideOutputSpec, ...] = (
    SlideOutputSpec(
        "ism_molecule_slide",
        "ISM/CSM detections slide",
        "Balanced PowerPoint slide of secure ISM/CSM molecules.",
        "astro_molecules_{view_choice}.pptx",
        write_molecule_slide,
        writer_kwargs=_balanced_slide_kwargs,
        layout_builder=build_molecule_slide_layout,
        report_label="ISM/CSM slide layout report",
        report_description="Text report for the generated ISM/CSM slide layout.",
        report_filename_template="astro_molecules_{view_choice}_layout.md",
    ),
    SlideOutputSpec(
        "ppd_detection_slide",
        "PPD detections slide",
        "PowerPoint slide of PPD molecules and isotopologues.",
        "ppd_molecules_{view_choice}.pptx",
        write_ppd_detection_slide,
        layout_builder=build_ppd_detection_slide_layout,
        report_label="PPD slide layout report",
        report_description="Text report for the generated PPD slide layout.",
        report_filename_template="ppd_molecules_{view_choice}_layout.md",
    ),
)


def figure_output_names() -> tuple[str, ...]:
    """Return stable names for standard figure outputs."""
    return tuple(spec.name for spec in FIGURE_OUTPUTS)


def table_output_names() -> tuple[str, ...]:
    """Return stable names for standard table output groups."""
    return tuple(spec.name for spec in TABLE_OUTPUTS)


def slide_output_names() -> tuple[str, ...]:
    """Return stable names for standard slide outputs."""
    return tuple(spec.name for spec in SLIDE_OUTPUTS)


__all__ = [
    "FIGURE_OUTPUTS",
    "SLIDE_OUTPUTS",
    "TABLE_OUTPUTS",
    "FigureOutputSpec",
    "OutputContext",
    "SlideOutputSpec",
    "TableOutputSpec",
    "figure_output_names",
    "slide_output_names",
    "table_output_names",
]
