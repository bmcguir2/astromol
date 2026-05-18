# astromol Architecture Specification

## Overview

astromol is a Python library for cataloging astronomical molecule detections.
Data is stored in JSON files, with bibliographic metadata loaded from a Zotero
BibTeX export. The code loads data, resolves cross-references, and provides
query helpers through `Database`.

## File Structure

```text
astromol/
├── data/
│   ├── molecules.json
│   ├── detections.json
│   ├── sources.json
│   ├── telescopes.json
│   └── references.bib
├── __init__.py
├── models.py
└── database.py
```

`test_load.py` is a local smoke script. The data files may be incomplete during
curation, but loaded cross-references must resolve for `Database()` to succeed.

`pyproject.toml` defines the installable package metadata, core dependencies,
optional dependency groups, package data, and test discovery configuration.

## Data Provenance

The current `molecules.json`, `detections.json`, and `sources.json` baseline was
promoted from the legacy conversion workflow in commit `67f37d9`. The final
audited preview state is preserved in commit `0187018`; use those commits for
the removed conversion scripts, staging reports, and legacy source data.

Molecule history records were imported from the legacy `census_version` and
`change_log` fields through a staged preview workflow. Legacy versions are
preserved as `history.events[*].legacy_version`, while
`history.introduced.census` records when a molecule entered astromol tracking
or census discussion. Confirmed census membership is tracked separately with
`history.accepted.census`, which distinguishes accepted entries from tentative
or disputed molecules discussed in an earlier census.

Molecule manuscript prose from the 2021 census LaTeX source was imported into
`molecules.json` through a staged preview workflow in commit `f2f4c4d`. That
commit preserves the one-off parser, source TeX, preview JSON, and review
report used for the audit trail; those staging artifacts are not part of the
working data model.

Historical census table verification was performed against the 2021 table
sources in commit `8e6b1dd` and the 2018 ISM/CSM, exoplanet, extragalactic,
and protoplanetary-disk tables in commit `d93b902`. Those commits preserve the
temporary table inputs and verification report used to validate accepted census
membership and detection-history boundaries. The table audit artifacts are not
part of the working data model.

Database-backed census-summary figure generation was migrated and verified in
commit `6c5c89f`. That checkpoint preserves the temporary 2021 arXiv reference
bundle, legacy molecule module, migrated `astromol.figures` implementation,
figure regression scripts, and migration documentation used for the audit
trail. The arXiv/legacy reference inputs were removed after that checkpoint and
are not part of the working data model.

OpenAI Codex is used as a coding assistant for software implementation,
refactoring, validation scripts, documentation drafts, and workflow support.
Scientific database content is not AI-generated. Molecule records, detection
records, source/telescope metadata, reference mappings, classifications, and
curation decisions are derived from the scientific literature, Zotero-managed
bibliography exports, legacy census materials, and human curator review. AI
assistance may be used to stage or transform records for review, but no staged
scientific data are accepted into production without human verification.

When a database value is calculated as part of this project, the supporting
notebook should be committed under `docs/calculations/`. Structured reference
fields still contain BibTeX citekeys, so they should cite the software, method,
basis-set, and/or literature sources that support the calculation. The local
notebook path should be recorded in a note or history summary for provenance.
Local `.ipynb_checkpoints/` directories are not tracked.

## Curation Workflow

Production data lives in the JSON files under `astromol/data/`. New records
should normally be staged first using copy/paste YAML templates from
`curation/templates/`, with working staging files kept under `curation/staging/`.
Run `python scripts/stage_records.py --staging <file.yaml>` to generate preview
JSON and a review report. Use `--apply` only after the preview is accepted.
The staging script automatically populates semantic history metadata for new
records, including `history.introduced.date`, `history.last_modified`, and an
initial `added` event dated with the staging run date. New molecule records and
secure detection records also default to a current-census `history.accepted`
block; set `history.accepted: null` for records that are tracked but not yet
accepted as confirmed.

Production inventory counts and accepted validation warnings are tracked in the
committed baseline `tests/baselines/production_data.json`. After an approved
data update is applied, run `python scripts/update_data_baseline.py` and inspect
the baseline diff before committing. Tests read this file so stale baselines
fail locally and in CI.

When the schema for `Molecule`, `Detection`, `Source`, or `Telescope` changes,
the corresponding curation template, `scripts/stage_records.py`, and this
specification must be updated in the same change.

## Census Views

Use `astromol.census.CensusView` as the shared filtering layer for manuscript
tables, figures, slides, scalar counts, and other census outputs. Output
generators should consume a view instead of filtering raw `Database` records
directly.

Two scope modes are supported:

- `CensusView.for_census(db, "2021")`: frozen/historical census boundary.
  Accepted secure records are selected with `history.accepted.census <= 2021`.
- `CensusView.current(db)`: live database boundary. Accepted secure records are
  selected regardless of census label.

Views exclude isotopologue records by default for molecule counts, context
tables, source/facility counts, and other standard census products. Pass
`include_isotopologues=True` to include isotopologue detections and molecules in
the same census/context query. This option is available for all census scopes
and all detection contexts.

Before the 2026 census cutoff, `CensusView.for_census(db, "2026")` and
`CensusView.current(db)` should return the same secure accepted detection sets
under the same isotopologue setting. After the 2026 census is frozen, they may
diverge as new post-2026 records are added to the live database.

Tentative and disputed detections are excluded by default. They can be included
for historical reproduction or review with `include_tentative=True` and
`include_disputed=True`, in which case the view uses record introduction
history rather than accepted history.

## Output Generation

Manuscript, figure, table, and slide generators should live in the package API
and consume `CensusView` instances. They must not duplicate census-boundary
filtering against raw JSON records.

`astromol.latex` contains the first migrated output helpers. Scalar manuscript
inputs are generated as filename-to-content mappings with `scalar_fragments`
and can be written with `write_scalar_fragments`. Generated scalar fragments
use the legacy `\endinput` convention so they can be included directly by
LaTeX manuscript sources. ISM/CSM molecule tables are generated with
`ism_table_fragments` and `write_ism_tables`. Use `layout="legacy"` for audited
2021-style reproduction and `layout="balanced"` for production manuscript
output. Membership comes from the view, table ordering uses first accepted
non-isotopologue ISM/CSM detection date within each atom-count category, and
balanced output enforces configurable row and column limits while preserving
explicit molecule labels through `\molref{label}{formula}` cells. Main ISM/CSM
molecule tables intentionally exclude isotopologue records. When balanced output
splits one atom-count category across multiple columns, the shared category
header is rendered with `\multicolumn` and centered over those columns.
Production ISM tables use a `13+ Atoms` terminal atom-count category before
PAHs and fullerenes, and should keep visual density constant or decreasing
across successive table fragments.

External-galaxy molecule tables are generated with `exgal_table_fragments` and
`write_exgal_table`. These tables include secure extragalactic detections and
introduced tentative extragalactic detections, mark tentative rows with a
dagger, number observation references by first appearance, and intentionally
exclude isotopologue records. Extragalactic atom-count table groups are derived
from non-empty bins rather than hard-coded to the current census inventory.

Protoplanetary-disk molecule tables are generated with `ppd_table_fragments`
and `write_ppd_table`. These tables intentionally include secure isotopologue
detections, exclude tentative/disputed detections, and order rows parent-first
so detected isotopologues are listed immediately after their parent molecule
where possible. PPD atom-count table groups are derived from non-empty bins.

Exoplanet-atmosphere molecule tables are generated with
`exoplanet_table_fragments` and `write_exoplanet_table`. The standard
exoplanet table excludes isotopologue records, excludes tentative/disputed
detections, and numbers observation references by first appearance. Pass
`include_isotopologues=True` only for an explicitly isotope-expanded exoplanet
output.

Interstellar-ice molecule tables are generated with `ice_table_fragments` and
`write_ice_table`. The standard ice table includes secure and tentative
detections, excludes isotopologue records, and marks tentative rows such as the
OCN- ice claim with a dagger. Pass `include_tentative=False` only for an
explicit secure-only output.

Detection-rate-by-atoms tables are generated with
`rate_by_atoms_table_fragments` and `write_rate_by_atoms_table`. They use
secure, non-isotopologue ISM/CSM first-detection years from the supplied
`CensusView`. The terminal normal-molecule size bin is a real `13+` category
that excludes PAHs and fullerenes, which are fit separately. Rates and `R^2`
values are computed from least-squares linear fits with NumPy.

Facility-count tables are generated with `facility_table_fragments` and
`write_facility_table`. They count observing-facility contributions from
secure, non-isotopologue ISM/CSM detections in the supplied `CensusView`.
Facilities with no detections in the view are omitted. Rows are sorted by
descending count with alphabetical tie-breaking and rendered as two side-by-side
facility/count column pairs.

Source-count tables are generated with `source_table_fragments` and
`write_source_table`. They count source contributions from secure,
non-isotopologue ISM/CSM detections in the supplied `CensusView`. Sources with
no detections in the view are omitted. Diffuse-cloud line-of-sight sources are
consolidated into a single `Diffuse Cloud` row for manuscript-table
reproduction, while closely related source-region grouping is inherited from
the curated `Source` records. Rows are sorted by descending count with
alphabetical tie-breaking and rendered as two side-by-side source/count column
pairs.

Figures are generated in `astromol.figures`. Figure modules should expose
testable data-builder functions separately from plotting/writing functions.
Manuscript figures should share the common figure-size, fixed axes rectangle,
axis-label, tick, and in-plot annotation defaults defined in `astromol.figures`
unless there is a specific scientific or layout reason to diverge. Figure
writers use the fixed canvas and axes rectangle rather than automatic tight
cropping when paired panels need matched plot areas. Current and future-facing
figures should use the package's color-blind-friendly palette by default;
legacy 2018/2021 reproduction views may preserve original colors when visual
parity is the goal. Plot axes and spines should be black by default, with a
full box around the plot area unless a specific figure design requires a
documented exception. Figure writers save raster output at 300 DPI, including
PNG files and any rasterized artists embedded in vector formats such as PDF.
The current figure API covers the database-generated census-summary figures
used in the 2021 manuscript body. Static appendix/conceptual graphics from the
2021 arXiv bundle are outside the database-output API unless they are rebuilt
as explicit future package products.
The default project display color is `ASTROMOL_BLUE` (`dodgerblue`). Use
`MIT_RED` (`#750014`) only for true emphasis/highlights, not as a general
secondary-series color. Use `NRAO_BLUE` (`#0A1589`) for dark outlines or
dark-blue elements, and `MIT_GRAY` (`#8B959E`) for secondary gray elements
only when gray is semantically appropriate. Black remains the standard color
for primary text, axes, ticks, and ordinary labels.
Multi-category scientific plots may still use the
color-blind-friendly palette when distinct categories need separate colors.
The standard cumulative-detections figure is built with
`cumulative_detection_data`, `plot_cumulative_detections`, and
`write_cumulative_detections_plot`; it uses secure, non-isotopologue ISM/CSM
first-detection years from the supplied `CensusView`. Cumulative-detection
trend annotations are census-aware: legacy 2018/2021 views retain facility
commissioning labels, while 2026/current views omit them. Trend annotations
and total-count labels are placed in axes-fraction coordinates with explicit
y-axis headroom so figure text remains stable as the database grows.
Cumulative-by-atom-count figures use the same secure, non-isotopologue ISM/CSM
first-detection records and group them into 2-12 atom traces, a true
non-PAH/non-fullerene `13+` trace, and separate PAH and fullerene traces.
`plot_cumulative_by_atoms` renders the legacy line view, while
`plot_stacked_cumulative_by_atoms` renders a composition-focused stacked view
using the same category colors and in-plot legend placement. The atom-count
palette is view-aware: 2018/2021 legacy views keep the original palette, while
2026/current and future views use the color-blind-friendly palette.
`rolling_rate_by_atoms_heatmap_data` and
`write_rolling_rate_by_atoms_heatmap` render the promoted rolling-rate
companion figure. The production heatmap uses secure, non-isotopologue ISM/CSM
first-detection records, the same 2-12, `13+`, fullerene, and PAH category
definitions, and a 10-year trailing detection-rate window by default. It uses
a compact single-column canvas with the same 10:8 aspect ratio as the other
cumulative atom-count figures, a `cividis` sequential palette, and short row
labels (`Fuller` denotes the fullerene class, not C60 specifically).
`detection_rate_by_atoms_data` provides the average-detection-rate bubble-plot
data. The rate is defined as the number of first detections in an atom-count
category divided by the number of years from that category's first detection
through the selected census/current end year. `write_detection_rate_by_atoms_plot`
preserves the single-view legacy reproduction path. The production 2026
version is `write_detection_rate_by_atoms_comparison_plot`, which overlays the
selected current/census view on the 2021 values. The comparison plot uses one
absolute marker-area scale for both datasets (`area = count * constant`) so
marker areas are directly comparable across census years. This deliberately
differs from the 2021 legacy plot, which normalized marker sizes within one
dataset. The production x-axis displays only 2-12 atom categories.
`facility_share_data` supplies facility contribution shares for both the
production bar chart and the legacy pie-grid reproduction. Facility shares use
secure, non-isotopologue ISM/CSM first-detection records. Each facility's
percentage is its number of first-detection contributions divided by all first
detections that occurred during that facility's operational window in the
selected view. The 2021 reproduction path preserves the published facility
set; current/future views select the top facilities from the curated database
by first-detection contribution count. The preferred production view is
`write_facility_share_bars_plot`, which labels each bar with `percent (n/N)`.
`write_facility_shares_plot` is retained for the legacy 3x3 pie-grid figure.
Active facilities use the project blue; inactive facilities use the legacy
pastel red (`#F87070`) because MIT Red is reserved for emphasis and is too
visually heavy for this role.
`scopes_by_year_data` and `write_scopes_by_year_plot` generate cumulative
first-detection contribution traces for facilities above a configurable
minimum contribution threshold. The standard view uses secure,
non-isotopologue ISM/CSM first detections and preserves the legacy fit cutoffs
for NRAO 36-ft, NRAO 140-ft, and Nobeyama 45-m when reporting per-year rates.
`plot_scopes_by_year(..., style="legacy")` preserves the legacy reproduction
style. `style="modern"` uses the 2026 production palette, emphasizes Yebes
40-m in MIT Red, and fades dotted post-last-contribution tails for NRAO 36-ft,
Nobeyama 45-m, and NRAO 140-ft.
`periodic_heatmap_data` and `write_periodic_heatmap` generate the
periodic-table elemental-composition heatmap from secure, non-isotopologue
ISM/CSM molecules by default.
`mass_by_wavelength_data` and `write_mass_by_wavelength_plot` generate the
molecular-mass KDE figure grouped by first-detection wavelength. The default
view uses secure, non-isotopologue ISM/CSM first detections and the legacy
bandwidth convention (`bw = 0.5`). The preferred 2026 manuscript view is
`write_mass_by_wavelength_boxplot`, which uses horizontal box/whisker summaries
with individual jittered detections. Production mass-by-wavelength views should
exclude fullerenes by default.
`molecules_by_wavelength_atoms_data` and
`write_molecules_by_wavelength_atoms_plot` preserve the legacy six-panel
KDE/histogram wavelength-by-atom-count figure. The preferred 2026 manuscript
view is `write_molecules_by_wavelength_atoms_bubble_heatmap`, which uses the
periodic-table yellow-to-red palette and encodes absolute counts of secure,
non-isotopologue, non-fullerene ISM/CSM first detections by wavelength category
and atom-count bin.
`du_histogram_data` and `write_du_histogram` preserve the legacy
degree-of-unsaturation histogram. The preferred 2026 manuscript view is
`write_du_bar_chart`, which renders exact DU values as discrete bars rather
than treating DU as a continuous distribution. DU views use secure,
non-isotopologue ISM/CSM molecules, exclude fullerenes by default, and only
include formulas containing H/D/N/C/O/S/F/Cl. The production bar chart excludes
negative-DU values from small protonated hydride/non-carbon species by default,
because those values are artifacts of applying the neutral/organic DU
expression outside its chemically meaningful domain.
`kappa_histogram_data` and `write_kappa_histogram` generate the Ray
asymmetry-parameter histogram for secure, non-isotopologue ISM/CSM molecules.
Linear rotors are treated as the perfectly prolate limiting case
`kappa = -1`; nonlinear molecules require all three rotational constants.
Molecules without a usable kappa are excluded from the plotted distribution.
The production histogram uses the standard manuscript aspect ratio, full boxed
axes, and a top guide arrow from prolate through asymmetric to oblate with
arrowheads anchored at `kappa = -1` and `kappa = +1`. Pass
`show_shape_guide=False` only when a bare histogram is needed.
`molecule_type_data` and `write_type_pie_chart` generate the concentric
molecule-type ring chart. The standard view uses secure, non-isotopologue
ISM/CSM molecules and counts each type category independently, so overlapping
classes such as neutral radicals, cyclic ions, PAHs, and fullerenes contribute
to every applicable ring. Production rendering sorts rings by descending
category count while keeping stable category colors; pass `order="legacy"` for
historical visual reproduction.
`source_type_data` and `write_source_pie_chart` generate the generalized
first-detection source-type ring chart. The standard view uses secure,
non-isotopologue ISM/CSM first detections and credits each molecule at most
once per generalized source category. Production rendering sorts rings by
descending category count while keeping stable category colors; pass
`order="legacy"` for historical visual reproduction.
`individual_source_data` and `write_individual_source_pie_chart` generate the
individual first-detection source ring chart. The standard view uses secure,
non-isotopologue ISM/CSM first detections. This figure counts listed source
contributions rather than unique molecules per named source: Sgr B2, TMC-1,
IRC+10216, Orion, and Other are credited independently, and the denominator is
the number of first-detected molecules. As a result, percentages can sum to
more than 100% when one first detection lists multiple sources. Production
rendering sorts rings by descending category count while keeping stable
category colors; pass `order="legacy"` for historical visual reproduction.
`molecule_type_by_source_type_data` and
`write_molecule_type_by_source_type` generate the four-panel pie-grid view of
molecule-type composition within generalized first-detection source categories.
The standard view uses secure, non-isotopologue ISM/CSM first detections. A
molecule is credited at most once per source category, while molecule-type
classes are counted independently within each source category to preserve the
overlap convention used for neutral radicals, cyclic ions, and related classes.
The preferred 2026 production view is
`write_molecule_type_by_source_enrichment_matrix`, which compares each
source/type fraction to the overall secure ISM/CSM molecule-type fraction.
Cells report the fractional enrichment factor and raw count; the color scale is
log-centered on `1x`, using gray for depletion and astromol blue for
enrichment.
`du_by_source_type_data` and `write_du_by_source_type` preserve the legacy
degree-of-unsaturation KDE comparison across generalized first-detection source
types. The standard view uses secure, non-isotopologue ISM/CSM first detections,
credits each molecule at most once per source category, excludes fullerenes by
default, and follows the legacy KDE bandwidth convention. Count labels are
placed automatically by default, with manual `count_label_overrides` available
for final manuscript tuning.
The preferred 2026 production view is `write_du_by_source_type_boxplot`, which
uses horizontal box/whisker summaries plus jittered exact DU values. The
production box/strip plot excludes negative-DU domain artifacts by default and
uses the same box, whisker, median, and point conventions as
`write_mass_by_wavelength_boxplot`. Production row order is dark cloud, carbon
star, SFR, and diffuse cloud. The DU box/strip x-axis includes a small left
padding so zero-valued detections do not sit directly on the y-axis spine.
`relative_du_by_source_type_data` and `write_relative_du_by_source_type`
preserve the legacy relative-DU KDE panel grouped by generalized source type.
Relative DU is defined as `du / maxdu`. The preferred 2026 production view is
`write_relative_du_by_source_type_boxplot`, which uses the same box/strip
conventions and row order as the absolute DU-by-source production view. The
relative-DU production view also excludes negative-DU domain artifacts by
default and pads the left x-axis edge so zero-valued points remain visible.
`mass_by_source_type_data` and `write_mass_by_source_type` preserve the legacy
molecular-mass KDE comparison across generalized first-detection source types.
The standard view uses secure, non-isotopologue ISM/CSM first detections,
credits each molecule at most once per source category, excludes fullerenes by
default, and follows the legacy KDE bandwidth convention. The preferred 2026
production view is `write_mass_by_source_type_boxplot`, which uses the same
horizontal box/whisker plus jittered-point conventions and source-row order as
the DU/source boxplots.
`wavelength_by_source_type_data` and `write_wavelength_by_source_type` generate
the wavelength-by-source-type pie grid. The standard view uses secure,
non-isotopologue ISM/CSM first detections, credits each molecule at most once
per generalized source category, and credits every wavelength listed for that
first detection. UV and visible detections are combined into a single `UV/Vis`
wedge in the displayed diffuse-cloud panel to preserve the legacy manuscript
convention.
The preferred 2026 production view is
`write_wavelength_by_source_type_stacked_bar`, which renders normalized
horizontal stacked bars from the same wavelength-credit counts and labels the
wavelength-credit denominator for each source category.

PowerPoint slide generation lives in `astromol.slides`. Slide generators should
follow the same selection/layout/rendering separation used by tables and
figures: `CensusView` supplies the scientifically selected records, a layout
plan records atom-bin membership and geometry, and `python-pptx` renders the
final binary slide. The first implemented slide product is the ISM/CSM
molecule summary slide. `selected_molecule_slide_entries` returns secure,
non-isotopologue ISM/CSM molecules by default, ordered by first accepted
detection represented in the view. `build_molecule_slide_layout` builds an
auditable layout object with group/column counts and warnings, while
`write_molecule_slide` writes the PowerPoint file. The `profile="legacy"`
layout preserves the hand-tuned 2021 `make_mols_slide` geometry for audit and
reproduction. The `profile="balanced"` layout derives column counts and group
spacing from the selected inventory and should be used for current-census
production review. The balanced layout uses an upper molecule band for 2-10
atom groups and a lower-right band for 11, 12, and 13+ atom groups. Groups 2-6
have a two-column minimum, and the shorter lower band has stricter
row-capacity behavior so large-molecule groups move into multiple columns
before they overflow. Slide generation must not silently clip or overlap
molecules when an inventory outgrows a layout; it should expose warnings or
fail clearly so the layout can be retuned.

Slide metadata separates software traceability from data currency. The
top-right credit block reports the installed `astromol` package version when
package metadata are available. When the package is run directly from a
development checkout, it reports `development` plus the current git short hash
and a dirty marker if the checkout has uncommitted changes. The molecule-count
footer reports the latest curated-record modification date from the selected
view and is not a package version.

`write_ppd_detection_slide` is the PPD-focused slide product. It uses the same
selection/layout/rendering machinery as the ISM/CSM molecule slide but selects
secure PPD records and includes detected isotopologues by default. Its default
`profile="compact"` renders only occupied atom-count groups, which avoids empty
ISM-style bins on sparse PPD inventories. The compact profile chooses a larger
readable molecule font, spreads occupied groups across the slide, and exposes
the same layout warnings if a future PPD inventory outgrows the available
geometry.

`astromol.cli` provides the first-pass public command-line interface for
standard outputs. The installed `astromol` command lists registry names and
generates individual standard figures, table groups, and slide decks from
`astromol.registry`. The `astromol outputs` subcommand delegates to the standard
bundle generator used by `astromol-generate-outputs`, preserving the existing
automation entry point while exposing a shorter public command.

## Packaging

`pyproject.toml` is the source of truth for package metadata and Python
dependencies. The base install includes the database loader, census views,
models, LaTeX/table helpers, figure generation, PowerPoint slide generation,
and bundled production data files. Output generation is a core package
capability, so its dependencies are installed by default. YAML curation
templates and staging scripts are contributor-facing source-checkout tools until
they are promoted into a packaged CLI. Optional extras are limited to
contributor-oriented tooling:

- `docs`: Sphinx/Read the Docs build dependencies
- `dev`: test, build, and package-check dependencies

The package version is the software/API version and is separate from census
boundaries and curated-record modification dates. `astromol.__version__`
resolves installed package metadata when available and falls back to the
development version string in source checkouts. The refactor development
version is `2026.0.0.dev0` so future releases sort after the legacy `astromol`
PyPI package (`2021.7`); this ordering constraint does not make package version
metadata the authoritative census-boundary record.

## Validation And Tests

`astromol.validation` contains semantic production-data checks that go beyond
basic `Database()` loading. It validates identifier namespaces, detection
relationship reciprocals, record-history ordering, spectroscopy numeric fields,
and `first` detection flag consistency. The installed command-line entry point
is:

```bash
astromol-validate
```

Source checkouts can run the same validator with:

```bash
python -m astromol.validation
```

Validation warnings are allowed for known unresolved curation follow-ups, such
as inherited legacy `*` dipole placeholders. Validation errors should block
release and data-update commits.

Regression tests live under `tests/` and are run with:

```bash
python -m pytest
```

The current pytest suite uses a wrapped regression-script harness to preserve
the audited migration checks while making them visible to pytest. This is a
bridge, not the final testing architecture. After the public API, documentation
examples, and CI workflow settle, convert the highest-value regression scripts
incrementally into conventional unit and integration tests with fixtures,
focused assertions, and slow-test markers for figure and slide generation.

## Data Model

### Ref

Bibliographic reference loaded from `references.bib`.

Fields:
- `bibcode` (str, required): BibTeX citation key using
  `LastName:Year:FirstPage`, e.g. `Swings:1937:483`. Add lowercase suffixes
  only when needed to disambiguate duplicate keys, e.g. `Cabezas:2021:L9a`.
- `author` (str, required): compact author display string
- `journal` (str, required): journal, booktitle, publisher, or entry type
- `year` (int, required)
- `month` (int or null)
- `day` (int or null)
- `volume` (str or null)
- `page` (str or null)
- `title` (str or null)
- `doi` (str or null)
- `note` (str or null): currently populated from BibTeX `annotation`

Computed:
- `sortdate`

### Telescope

Observing facility or instrument.

Required fields:
- `name`
- `nick`
- `shortname`
- `type`
- `wavelength` (list of values from `WAVELENGTHS`)

Optional fields:
- `diameter`
- `latitude`
- `longitude`
- `built`
- `decommissioned`
- `note`
- `latex_name`
- `history`: optional `RecordHistory`

Computed:
- `active`

### Source

Astronomical object where molecules have been detected.

Required fields:
- `name`
- `nick`
- `type` (one of `SOURCE_TYPES`)

Optional fields:
- `ra`
- `dec`
- `simbad_url`
- `latex_name`
- `note`
- `history`: optional `RecordHistory`

Current `SOURCE_TYPES`:
- `SFR`
- `Dark Cloud`
- `Carbon Star`
- `Oxygen Star`
- `AGB Star`
- `Planetary Nebula`
- `Protostellar`
- `Diffuse Cloud`
- `HII Region`
- `PDR`
- `Shock`
- `Supernova Remnant`
- `Sgr A`
- `External Galaxy`
- `YSO`
- `Exoplanet`
- `Other`

### RotationalConstants

Rotational constants in MHz.

Fields:
- `A`
- `B`
- `C`
- `refs` (list): BibTeX keys in JSON, resolved to `Ref` objects during loading
- `note`

Compatibility:
- `ref` returns the first item in `refs`, or `None`.

### DipoleMoment

Dipole moment components in Debye.

Fields:
- `a`
- `b`
- `c`
- `refs` (list): BibTeX keys in JSON, resolved to `Ref` objects during loading
- `note`

Computed:
- `total`

Compatibility:
- `ref` returns the first item in `refs`, or `None`.

### RecordHistory

Semantic update history exposed through database records. This is not intended
to replace git as the fine-grained audit trail. It records user-facing database
events that explain when and why a record was added, updated, or corrected.

History fields:
- `introduced`: dict with optional `date`, `census`, and `context` keys. This
  records when a molecule or record entered astromol tracking or census
  discussion.
- `accepted`: dict, or `null`, with optional `date`, `census`, and `context`
  keys. For molecules, this records when the molecule first became an accepted
  confirmed census entry. For detections, this records when that detection
  context first became part of an accepted census table or census update. Use
  `null` for tentative or disputed records.
- `last_modified`: ISO date string for the last semantic database change
- `events`: list of `HistoryEvent`

`HistoryEvent` fields:
- `kind`: one of `added`, `updated`, or `corrected`
- `summary`: concise human-readable description of the event
- `date`: optional ISO date string
- `census`: optional census release string, such as `2021` or `2026`
- `fields`: optional list of affected field names
- `legacy_version`: optional legacy astromol version identifier preserved from
  the pre-refactor database

Molecules accepted in a general census table should be selected with molecule
`history.accepted.census`, not `history.introduced.census` and not detection
year alone. Context-specific tables, such as ice, PPD, extragalactic, and
exoplanet tables, should be selected with detection
`history.accepted.census`. This matters for records discussed as tentative or
disputed in an earlier census, confirmed later, or published before a census
cutoff but not accepted into that census table. Records changed for a given
census can be selected from `history.events[*].census`, independently of
whether they were newly introduced or newly accepted.

### Molecule

Intrinsic chemical species metadata.

Required fields:
- `name`
- `formula`: canonical internal chemical formula. This must be parseable by
  `molmass`; isotopologues use explicit bracket notation such as `H2[13C]O`.

Display and notes:
- `table_formula`: readable formula body for manuscript tables. This may use
  mhchem isotope placement such as `H2^{13}CO`, but it does not include the
  surrounding `\ce{}` wrapper; LaTeX generation adds that wrapper.
  Defaults to `formula`.
- `label` (required by `__post_init__`)
- `note`

Optional identifiers:
- `iupac_name`
- `selfies`
- `synonyms`
- `smiles`
- `canonical_smiles`
- `inchi`
- `inchikey`

Curated flags and structure metadata:
- `radical_override`
- `fullerene`
- `pah`
- `n_rings`
- `cyclic`
- `tags`: dictionary of manually curated analysis annotations. Values are lists
  of strings. Current category:
  - `functional_groups`: manually assigned functional-group or structural-motif
    tags used for analysis. These are curated annotations, not computed
    structural assertions.

Functional-group tag definitions:
- `CN`: Cyano/nitrile-family C-N motif. Includes nitriles, metal cyanides,
  cyanopolyynes, and related species manually classified as containing a
  cyano-like C-N unit.
- `NC`: Isocyano/isocyanide-family N-C motif, where the connectivity is the
  structural isomer of a cyano/nitrile-like group.
- `CH3`: Methyl group or methyl substituent.
- `=O`: Carbonyl/oxygen multiple-bond motif, broadly including aldehydes,
  ketones, acids, esters, amides, ureas, ketenes, isocyanates, and
  O-terminated carbon-chain analogues.
- `=S`: Thiocarbonyl/sulfur multiple-bond motif, broadly the sulfur analogue
  of `=O`, including thioaldehydes, thioketones, thioketenes,
  isothiocyanates, and S-terminated carbon-chain analogues.
- `OH`: Hydroxyl group, including alcohols, hydroxy-substituted species,
  carboxylic acids, and metal hydroxides.
- `SH`: Thiol/mercapto group, including sulfur analogues of hydroxyl-bearing
  species.
- `NH2`: Amino group or primary amide/amine-like `NH2` substituent.
- `NH`: Imine/imino/secondary amine-like `NH` motif, manually assigned where
  `NH` is structurally meaningful rather than just present in the formula.

Spectroscopy:
- `rotcon` (`RotationalConstants` or null)
- `dipole` (`DipoleMoment` or null)

References and relationships:
- `refs`: role-keyed dict using roles from `MOLECULE_REF_ROLES`
  - `lab` includes laboratory characterization references for the molecule.
    Legacy context-specific lab references, such as ice laboratory spectra, are
    normalized here rather than attached to detections.
  - `computation` includes computational-method, software, basis-set, or
    literature references supporting calculated molecular properties. If a value
    is computed as part of this project, record the calculation notebook path in
    the relevant property note or history entry because local notebooks are not
    BibTeX citekeys.
- `isotopologue_of`: parent molecule label or null

LaTeX fields:
- `latex_section_override`: optional custom section heading. This should be
  null for normal generated section headings.
- `latex_body`: curated prose for the molecule's manuscript section.

History:
- `history`: optional `RecordHistory`

Computed from formula or spectroscopy:
- `atom_counts`
- `atoms`
- `isotope_counts`
- `mass`
- `average_mass`
- `nominal_mass`
- `natoms`
- `charge`
- `cation`
- `anion`
- `neutral`
- `nelectrons`
- `odd_electron`
- `radical`
- `du`
- `maxdu`
- `is_linear`
- `kappa`

When a formula uses an isotope label that is valid for the database but
unsupported by `molmass`, the model uses `ISOTOPE_MASS_OVERRIDES` for hardcoded
relative atomic masses. Current override:
- `26Al`: 25.98689188 u

Loader compatibility:
- Legacy stored/computed fields such as `mass`, `natoms`, `charge`, `du`,
  `maxdu`, `is_linear`, and `kappa` are ignored before constructing
  `Molecule`.
- Legacy `radical` is mapped to `radical_override`.
- Legacy flat spectroscopy keys `Acon/Bcon/Ccon` and `mua/mub/muc` are mapped
  into `rotcon` and `dipole`.

### Detection

Specific molecule detection event.

Required fields:
- `id`: stable detection identifier
- `molecule`: molecule label, resolved to a `Molecule`
- `sources`: source nicks, resolved to `Source` objects
- `telescopes`: telescope nicks, resolved to `Telescope` objects
- `wavelengths`: values from `WAVELENGTHS`
- `year`
- `type`: value from `DETECTION_TYPES`

Optional fields:
- `note`
- `status`: value from `DETECTION_STATUSES`; defaults to `secure`
- `status_note`: concise explanation for tentative or disputed detections
- `first`
- `month`
- `day`
- `refs`: role-keyed dict using roles from `DETECTION_REF_ROLES`
- `confirms`: detection IDs this record confirms
- `confirmed_by`: detection IDs that confirm this record
- `disputes`: detection IDs this record disputes
- `disputed_by`: detection IDs that dispute this record
- `supersedes`: detection IDs this record supersedes
- `superseded_by`: detection IDs that supersede this record
- `latex_text`
- `history`: optional `RecordHistory`

Computed:
- `sortdate`

Detection IDs are stable database identifiers. The default format is:

```text
det:<molecule-label-without-mol-prefix>:<context>:<year>[:qualifier]
```

Examples:
- `det:CH:ism-csm:1937`
- `det:CH3CH2CCH:ism-csm:2021`
- `det:CH3CH2CCH:ism-csm:2024`
- `det:OCN-:ice:2005`

Use an additional qualifier only when molecule, context, and year would not be
unique. Detection IDs should not include mutable status terms such as
`tentative`, because a detection claim may later be confirmed without becoming a
different database record.

Detection status is independent of detection type. `type` records the
astrophysical context, such as `ISM/CSM`, `ice`, or `exgal`; `status` records
whether the claim is `secure`, `tentative`, or `disputed`. Molecule inclusion
does not imply a secure astronomical detection. Tentative and disputed records
should normally keep `first: false` unless a deliberate curation decision says
otherwise.

`first: true` means the first accepted/secure detection in that context, not
necessarily the first chronological claim. Earlier tentative or disputed claims
should be represented as separate detection records and linked to later secure
detections with `confirmed_by`/`confirms` where appropriate.

Relationship fields store detection IDs, not BibTeX keys. Use `refs` for papers
and relationship fields for links between detection records. For example, an
earlier tentative butyne detection can use
`confirmed_by: [det:CH3CH2CCH:ism-csm:2024]`, while the confirming paper remains
under `refs.confirmation` or the secure detection's `refs.observation`.

Current `DETECTION_TYPES`:
- `ISM/CSM`
- `isotopologue`
- `ice`
- `exgal`
- `ppd`
- `exo`
- `comet`

Current `DETECTION_STATUSES`:
- `secure`
- `tentative`
- `disputed`

Current `DETECTION_REF_ROLES`:
- `observation`: paper reporting the detection claim
- `confirmation`: paper confirming an earlier tentative detection
- `dispute`: paper challenging a detection claim
- `correction`: paper correcting detection metadata or interpretation

Current detection relationship fields:
- `confirms`
- `confirmed_by`
- `disputes`
- `disputed_by`
- `supersedes`
- `superseded_by`

## Reference Loading

`Database._load_refs()` reads `astromol/data/references.bib` using
`bibtexparser` and stores `Ref` objects in `Database.refs`, keyed by the BibTeX
citation key.

JSON references should use those BibTeX keys. The loader resolves:
- molecule-level `refs`
- detection-level `refs`
- nested spectroscopy `ref` or `refs`

Nested spectroscopy supports both current `refs: [...]` and legacy `ref`.

## Database

`Database()` loads data in this order:
1. references
2. telescopes
3. sources
4. molecules
5. detections

Primary stores:
- `refs`: dict keyed by BibTeX citation key
- `telescopes`: dict keyed by telescope nick
- `sources`: dict keyed by source nick
- `molecules`: dict keyed by unique molecule label
- `molecules_by_formula`: dict keyed by chemical formula, with list values
- `molecules_by_name`: dict keyed by molecule name, with list values
- `detections`: list of `Detection`

Molecule labels are the database identity keys. Labels should use the `mol:`
namespace, e.g. `mol:CH` or `mol:l-C3H`. Chemical formulae are data fields and
are not required to be unique because structural isomers can share a formula.
Detection `molecule` values and `isotopologue_of` values refer to molecule
labels, not formulae.

Accessors:
- `get_ref`
- `get_telescope`
- `get_source`
- `get_molecule`
- `get_molecules_by_formula`
- `get_molecules_by_name`

`Database.__repr__()` prints a compact count summary of loaded molecules,
detections, sources, telescopes, and references.

## Validation Expectations

For a complete dataset:
- all JSON files parse
- required dataclass fields are present
- all source types, detection types, and wavelengths are allowed
- molecule, source, telescope, and reference cross-references resolve
- molecule labels, source nicks, telescope nicks, and BibTeX keys are unique
- all `refs` entries use valid role names
