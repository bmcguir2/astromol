# Census Output Recipes

This page is the practical recipe book for regenerating census products from
the database. The same pattern is used throughout:

1. Load the production database.
2. Choose a `CensusView`.
3. Build the table/figure/slide data.
4. Write the output to a directory you control.

```python
from pathlib import Path

from astromol.census import CensusView
from astromol.database import Database

db = Database()
view = CensusView.for_census(db, "2026")
out = Path("census_outputs")
out.mkdir(exist_ok=True)
```

Use `CensusView.for_census(db, "2021")` for 2021 reproduction,
`CensusView.for_census(db, "2026")` for the developing 2026 census, and
`CensusView.current(db)` for the live database after a census boundary has
passed.

Until the 2026 cutoff is frozen, the 2026 and current views are expected to
select the same records.

## Latest Standard Output Bundle

The project publishes automatically generated standard figures and slides at
<https://bmcguir2.github.io/astromol/>. The same bundle can be generated locally
with:

```bash
astromol-generate-outputs --output-dir build/astromol_outputs --view current --formats png pdf
```

The command writes figures, PowerPoint slides, LaTeX table fragments, a static
`index.html`, and `astromol_latest_outputs.zip`. The landing page highlights the
bundle and the most commonly requested standard downloads first, then presents
public-facing figure download cards. LaTeX fragments and layout reports remain
available in the complete bundle.

For selective generation, use the public `astromol` command with the same
registry names:

```bash
astromol list figures
astromol figure cumulative_detections --view current --output cumulative_detections.pdf
astromol table ism_tables --view 2026 --output-dir build/tables
astromol slide ism_molecule_slide --view current --output-dir build/slides --report
```

The `outputs` subcommand is equivalent to the full-bundle workflow:

```bash
astromol outputs --output-dir build/astromol_outputs --view current --formats png pdf
```

## Standard View Options

Most table, figure, and slide data builders use the same filtering options:

| Option | Default | Meaning |
| --- | --- | --- |
| `detection_type` | `"ISM/CSM"` | Context to analyze. Other values include `"exgal"`, `"ppd"`, `"ice"`, and `"exo"`. |
| `include_tentative` | `False` | Include tentative detections by introduced history rather than accepted history. |
| `include_disputed` | `False` | Include disputed detections by introduced history. |
| `include_isotopologues` | `False` | Include molecule records marked as isotopologues. |
| `include_fullerenes` | varies | Include fullerenes in analyses where their extreme values can dominate the display. |

The manuscript defaults are intentionally conservative: secure detections only,
isotopologue-free for most ISM/CSM counts and tables, and fullerene exclusion
where the output is designed to compare ordinary molecular distributions.

## Recreate The Manuscript Tables

The LaTeX helpers live in `astromol.latex`. Writer functions create `.tex`
fragments and return the generated text as a dictionary.

```python
from astromol.latex import (
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

table_dir = out / "tables"
table_dir.mkdir(exist_ok=True)

write_scalar_fragments(view, table_dir)
write_ism_tables(view, table_dir, layout="balanced")
write_exgal_table(view, table_dir)
write_ppd_table(view, table_dir)
write_exoplanet_table(view, table_dir)
write_ice_table(view, table_dir)
write_rate_by_atoms_table(view, table_dir)
write_facility_table(view, table_dir)
write_source_table(view, table_dir)
```

### Table Defaults

| Helper | Default output | Main options |
| --- | --- | --- |
| `write_scalar_fragments` | Individual scalar `.tex` inputs such as molecule counts, source counts, and detection rates. | None beyond the view. |
| `write_ism_tables` | Balanced ISM/CSM molecule tables, isotopologue-free, ordered by first detection. | `layout="balanced"` or `"legacy"`, `max_rows`, `max_columns`. |
| `write_exgal_table` | External-galaxy molecule table, isotopologue-free; tentative detections are included and marked. | Filename/output directory. |
| `write_ppd_table` | PPD molecule table including isotopologues by default. | Filename/output directory. |
| `write_exoplanet_table` | Exoplanet-atmosphere table, secure detections only. | `include_isotopologues`. |
| `write_ice_table` | Ice table, tentative detections included and marked by default. | `include_tentative`, `include_isotopologues`. |
| `write_rate_by_atoms_table` | Detection-rate fit table by atom-count category. | Filename/output directory. |
| `write_facility_table` | Facility contribution count table. | Filename/output directory. |
| `write_source_table` | Source contribution count table with diffuse-cloud/LOS detections grouped. | Filename/output directory. |

To reproduce the historical 2021 ISM table geometry:

```python
view_2021 = CensusView.for_census(db, "2021")
write_ism_tables(view_2021, table_dir, layout="legacy", basename="ism_table_2021")
```

## Recreate The Manuscript Figures

Figure helpers live in `astromol.figures`. Most figures have one data-builder
function and one writer function. Build the data first, then pass it to the
writer.

```python
from astromol.figures import (
    cumulative_detection_data,
    write_cumulative_detections_plot,
)

data = cumulative_detection_data(view)
write_cumulative_detections_plot(data, out / "cumulative_detections.pdf")
```

The writer infers the output format from the file suffix. Use `.pdf` for
manuscript graphics and `.png` for quick previews or slides.

For standard production outputs, `astromol.registry.FIGURE_OUTPUTS` is the
canonical list. This is the same registry used by `astromol-generate-outputs`:

```python
from astromol.registry import FIGURE_OUTPUTS, OutputContext

context = OutputContext(
    view=view,
    view_choice="2026",
    baseline_view=CensusView.for_census(db, "2021"),
)

for spec in FIGURE_OUTPUTS:
    print(spec.name, "-", spec.label)
    spec.write(context, figure_dir / f"{spec.name}.pdf")
```

### Figure Recipe List

```python
from astromol.figures import (
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
    write_du_by_source_type,
    write_du_by_source_type_boxplot,
    write_du_histogram,
    write_facility_share_bars_plot,
    write_facility_shares_plot,
    write_individual_source_pie_chart,
    write_kappa_histogram,
    write_mass_by_source_type,
    write_mass_by_source_type_boxplot,
    write_mass_by_wavelength_boxplot,
    write_mass_by_wavelength_plot,
    write_molecule_type_by_source_enrichment_matrix,
    write_molecule_type_by_source_type,
    write_molecules_by_wavelength_atoms_bubble_heatmap,
    write_molecules_by_wavelength_atoms_plot,
    write_periodic_heatmap,
    write_relative_du_by_source_type,
    write_relative_du_by_source_type_boxplot,
    write_rolling_rate_by_atoms_heatmap,
    write_scopes_by_year_plot,
    write_source_pie_chart,
    write_stacked_cumulative_by_atoms_plot,
    write_type_pie_chart,
    write_wavelength_by_source_type,
    write_wavelength_by_source_type_stacked_bar,
)

figure_dir = out / "figures"
figure_dir.mkdir(exist_ok=True)

write_cumulative_detections_plot(
    cumulative_detection_data(view),
    figure_dir / "cumulative_detections.pdf",
)

by_atoms = cumulative_by_atoms_data(view)
write_cumulative_by_atoms_plot(by_atoms, figure_dir / "cumulative_by_atoms.pdf")
write_stacked_cumulative_by_atoms_plot(
    by_atoms,
    figure_dir / "cumulative_by_atoms_stacked.pdf",
)

write_rolling_rate_by_atoms_heatmap(
    rolling_rate_by_atoms_heatmap_data(view),
    figure_dir / "rolling_rate_by_atoms.pdf",
)

write_periodic_heatmap(
    periodic_heatmap_data(view),
    figure_dir / "periodic_heatmap.pdf",
)

du_data = du_histogram_data(view)
write_du_histogram(du_data, figure_dir / "du_histogram_legacy.pdf")
write_du_bar_chart(du_data, figure_dir / "du_bar_chart.pdf")

write_kappa_histogram(
    kappa_histogram_data(view),
    figure_dir / "kappa_distribution.pdf",
)

write_type_pie_chart(
    molecule_type_data(view),
    figure_dir / "molecule_type_pie.pdf",
)
write_source_pie_chart(
    source_type_data(view),
    figure_dir / "source_type_pie.pdf",
)
write_individual_source_pie_chart(
    individual_source_data(view),
    figure_dir / "individual_source_pie.pdf",
)

type_source = molecule_type_by_source_type_data(view)
write_molecule_type_by_source_type(
    type_source,
    figure_dir / "molecule_type_by_source_legacy.pdf",
)
write_molecule_type_by_source_enrichment_matrix(
    type_source,
    figure_dir / "molecule_type_source_enrichment.pdf",
)

du_source = du_by_source_type_data(view)
write_du_by_source_type(du_source, figure_dir / "du_by_source_legacy.pdf")
write_du_by_source_type_boxplot(
    du_source,
    figure_dir / "du_by_source_boxplot.pdf",
)

relative_du_source = relative_du_by_source_type_data(view)
write_relative_du_by_source_type(
    relative_du_source,
    figure_dir / "relative_du_by_source_legacy.pdf",
)
write_relative_du_by_source_type_boxplot(
    relative_du_source,
    figure_dir / "relative_du_by_source_boxplot.pdf",
)

mass_source = mass_by_source_type_data(view)
write_mass_by_source_type(
    mass_source,
    figure_dir / "mass_by_source_legacy.pdf",
)
write_mass_by_source_type_boxplot(
    mass_source,
    figure_dir / "mass_by_source_boxplot.pdf",
)

wavelength_source = wavelength_by_source_type_data(view)
write_wavelength_by_source_type(
    wavelength_source,
    figure_dir / "wavelength_by_source_legacy.pdf",
)
write_wavelength_by_source_type_stacked_bar(
    wavelength_source,
    figure_dir / "wavelength_by_source_stacked.pdf",
)

mass_wave = mass_by_wavelength_data(view)
write_mass_by_wavelength_plot(
    mass_wave,
    figure_dir / "mass_by_wavelength_legacy.pdf",
)
write_mass_by_wavelength_boxplot(
    mass_wave,
    figure_dir / "mass_by_wavelength_boxplot.pdf",
)

wave_atoms = molecules_by_wavelength_atoms_data(view)
write_molecules_by_wavelength_atoms_plot(
    wave_atoms,
    figure_dir / "wavelength_atoms_legacy.pdf",
)
write_molecules_by_wavelength_atoms_bubble_heatmap(
    wave_atoms,
    figure_dir / "wavelength_atoms_bubble_heatmap.pdf",
)

rate_current = detection_rate_by_atoms_data(view)
write_detection_rate_by_atoms_plot(
    rate_current,
    figure_dir / "detection_rate_by_atoms.pdf",
)

view_2021 = CensusView.for_census(db, "2021")
rate_2021 = detection_rate_by_atoms_data(view_2021)
write_detection_rate_by_atoms_comparison_plot(
    rate_current,
    rate_2021,
    figure_dir / "detection_rate_by_atoms_comparison.pdf",
    current_label="2026",
    baseline_label="2021",
)

facility_data = facility_share_data(view)
write_facility_shares_plot(
    facility_data,
    figure_dir / "facility_shares_legacy.pdf",
)
write_facility_share_bars_plot(
    facility_data,
    figure_dir / "facility_share_bars.pdf",
)

write_scopes_by_year_plot(
    scopes_by_year_data(view),
    figure_dir / "scopes_by_year.pdf",
    style="modern",
)
```

### Figure Defaults And Common Options

Standard figure writers save PNG output, and rasterized artists embedded in PDF
output, at 300 DPI.

| Figure | Default interpretation | Useful options |
| --- | --- | --- |
| Cumulative detections | Secure, non-isotopologue ISM/CSM first detections through the view boundary. | `detection_type`, `start_year`, `end_year`, tentative/disputed/isotopologue flags. |
| Cumulative by atoms | Atom-count traces: 2-12, 13+, PAHs, fullerenes. | `series_specs`, `include_isotopologues`. |
| Rolling rate heatmap | Trailing rolling detections/year by atom-count category. | `window`, `vmax`, `series_specs`. |
| Periodic heatmap | Number of selected molecules containing each element. | `detection_type`, tentative/disputed/isotopologue flags. |
| DU distribution | Degree of unsaturation for compatible formulas. | `include_fullerenes`, `include_isotopologues`. |
| Kappa distribution | Ray asymmetry parameter for molecules with usable rotational constants. | `detection_type`, tentative/disputed/isotopologue flags. |
| Molecule/source pie charts | Legacy categorical summaries. | `order`, context flags. |
| Molecule type/source enrichment | Relative over- or under-representation of molecule types by first-detection source category. | `source_label_overrides`. |
| DU/source, relative-DU/source, mass/source | Source-category comparisons. Legacy KDE and production box/strip versions are both available. | `include_fullerenes`, `source_label_overrides`. |
| Wavelength/source | Wavelength contributions by first-detection source category. | `include_fullerenes`, `source_label_overrides`. |
| Mass/wavelength | Mass distributions by first-detection wavelength. | `label_mode`, `label_overrides`; boxplot writer for production view. |
| Wavelength/atom-count | Molecule atom counts by detection wavelength. | Bubble heatmap writer for production view. |
| Detection rate by atoms | Average detections/year by atom-count category. | `categories`, `end_year`. |
| Facility shares | Facility contribution share during each facility window. | `top_n`, `facility_nicks`, `use_legacy_2021_selection`. |
| Scopes by year | Cumulative facility contribution traces. | `min_detections`, `cutoffs`, `colors_by_shortname`, `style`. |

## Recreate The Slides

PowerPoint helpers live in `astromol.slides`.

```python
from astromol.slides import (
    build_molecule_slide_layout,
    write_molecule_slide,
    write_molecule_slide_report,
    write_ppd_detection_slide,
)

slide_dir = out / "slides"
slide_dir.mkdir(exist_ok=True)

write_molecule_slide(
    view,
    slide_dir / "astro_molecules_2026.pptx",
    profile="balanced",
)

write_ppd_detection_slide(
    view,
    slide_dir / "ppd_molecules_2026.pptx",
)
```

### Slide Defaults

| Helper | Default output | Main options |
| --- | --- | --- |
| `write_molecule_slide` | ISM/CSM molecule slide. The default profile is `legacy`, but the current-census production profile is `balanced`. | `profile`, `detection_type`, tentative/disputed/isotopologue flags, `title`, `last_updated`. |
| `write_ppd_detection_slide` | PPD molecule/isotopologue slide. | Includes isotopologues by default; `profile="compact"`. |
| `build_molecule_slide_layout` | Layout plan without writing a `.pptx`. | Same selection options as `write_molecule_slide`. |
| `write_molecule_slide_report` | Markdown report of slide group counts, font size, and layout warnings. | Pass a layout from `build_molecule_slide_layout`. |

Generate a layout report before writing a slide:

```python
layout = build_molecule_slide_layout(view, profile="balanced")
write_molecule_slide_report(layout, slide_dir / "astro_molecules_layout.md")

if layout.warnings:
    for warning in layout.warnings:
        print(warning)
```

Historical 2021 reproduction:

```python
write_molecule_slide(
    view_2021,
    slide_dir / "astro_molecules_2021_legacy.pptx",
    profile="legacy",
)
```

## Custom Database Views

Use `CensusView.filtered(...)` when a standard census/current view has the
right historical boundary, but you want to run an output on a special subset.
The filters are applied after the normal census/status/isotopologue selection
rules.

### Carbon-Bearing ISM/CSM Molecules

```python
carbon_view = view.filtered(
    molecule_filter=lambda molecule: molecule.atom_counts.get("C", 0) > 0
)

write_periodic_heatmap(
    periodic_heatmap_data(carbon_view),
    figure_dir / "carbon_only_periodic_heatmap.pdf",
)
```

### Millimeter First Detections

```python
mm_view = view.filtered(
    detection_filter=lambda detection: "mm" in detection.wavelengths
)

write_cumulative_detections_plot(
    cumulative_detection_data(mm_view),
    figure_dir / "mm_cumulative_detections.pdf",
)
```

### Detections In One Source Type

```python
dark_cloud_view = view.filtered(
    detection_filter=lambda detection: any(
        source.type == "Dark Cloud"
        for source in detection.sources
    )
)

write_du_bar_chart(
    du_histogram_data(dark_cloud_view),
    figure_dir / "dark_cloud_du_bar_chart.pdf",
)
```

### Combining Filters

```python
carbon_dark_cloud_view = view.filtered(
    molecule_filter=lambda molecule: molecule.atom_counts.get("C", 0) > 0,
    detection_filter=lambda detection: any(
        source.type == "Dark Cloud"
        for source in detection.sources
    ),
)
```

Use custom filtered views for exploratory or supplemental analysis. For formal
census reproduction, use the standard census views and document any option
changes explicitly.
