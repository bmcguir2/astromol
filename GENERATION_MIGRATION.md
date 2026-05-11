# Census Output Generation Migration

This log tracks modernization of the legacy census-output functions from
`astromol/functions_legacy.py` into the refactored database-backed package API.
It is a working audit artifact for the migration; production scientific data
remain the JSON files under `astromol/data/`.

## Reference Inputs

- `astromol/functions_legacy.py`: legacy 2021 output-generation code.
- Commit `6c5c89f`: figure-migration audit checkpoint preserving the temporary
  2021 arXiv source/assets under `astromol/arxiv/` and
  `astromol/data/molecules_legacy.py`.

The legacy module is treated as a behavioral specification, not as runnable
production code. It imports removed legacy modules (`astromol.molecules`,
`astromol.sources`, `astromol.telescopes`) and contains many embedded LaTeX
strings that emit Python syntax warnings under modern parsing.

The temporary arXiv bundle and legacy molecule module were removed after the
audit checkpoint; recover them from commit `6c5c89f` if a future verification
pass needs the exact staged reference inputs.

## Decisions

- Manuscript-facing reminders for the 2026 paper are collected in
  `MANUSCRIPT_NOTES_2026.md`; use that file rather than scattering writing
  notes through migration or schema documentation.
- Output-generation code should live inside the installable `astromol` API.
- Reusable census-selection logic will be centralized before porting output
  functions. Tables, figures, and slides must consume the same view/filtering
  layer rather than independently filtering raw records.
- Users should eventually be able to generate an up-to-date complete compiled
  PDF, figures, tables, and the summary slide from the current database between
  formal publication releases.
- During the documentation phase, Read the Docs pages, example notebooks, and
  any Google Colab materials need clear, copy-pasteable commands for generating
  every supported figure, table, slide, and manuscript artifact from the
  database. This should include both package/API examples and command-line
  entry points once the CLI is defined.
- The PowerPoint molecule slide is a first-class product and will be ported.
- A PowerPoint slide summarizing PPD detections should also be added as a
  first-class product.
- The migration will proceed one function group at a time, with parity checks
  against the 2021 behavior where possible.
- Modern figure APIs should keep scientifically useful flexibility, especially
  for census views and data slices, but should not preserve every legacy
  cosmetic option. Prefer readable, explicit plotting code with named data
  preparation functions over compact code or highly parameterized visual
  wrappers.
- Every figure refinement pass should explicitly assess color-blind
  accessibility, small-format readability, line/marker distinguishability, and
  whether the figure is carrying too many visual categories for the intended
  manuscript use.
- Figure styling now has named project colors in `astromol.figures`: default
  display blue (`ASTROMOL_BLUE`, matplotlib `dodgerblue`), MIT Red
  (`#750014`) for true emphasis/highlights only, NRAO Blue (`#0A1589`) for
  dark outlines/elements, and MIT Gray (`#8B959E`) for secondary gray
  elements. Use the established color-blind-friendly palette instead when a
  figure needs many distinct category colors.

## Proposed Package Structure

- `astromol/census.py`: `CensusView` and shared census/context filters.
- `astromol/analysis.py`: reusable statistics and aggregation helpers.
- `astromol/latex.py`: scalar and table LaTeX generators.
- `astromol/figures.py`: figure data builders and plotting functions.
- `astromol/slides.py`: PowerPoint/slide generation.
- `astromol/manuscript.py`: high-level manuscript-output orchestration.

The exact names can change during implementation, but the ownership boundaries
should remain: filtering in one place, aggregations reusable, output renderers
thin.

## CensusView Target

The first implementation milestone is a shared view API, conceptually:

```python
from astromol.database import Database
from astromol.census import CensusView

db = Database()
view = CensusView(db, census="2021")

view.ism_detections()
view.ism_molecules()
view.ism_molecules(include_isotopologues=True)
view.ppd_detections()
view.exgal_detections(include_tentative=True)
view.exoplanet_detections()
view.source_counts()
view.facility_counts()
```

The core rules should use `history.accepted.census` for accepted membership,
not detection year alone. Tentative and disputed records remain queryable but
are excluded from accepted tables unless explicitly requested.

## Migration Phases

| Phase | Scope | Status |
| --- | --- | --- |
| 0 | Inventory legacy functions and classify target destinations. | Complete |
| 1 | Implement and verify `CensusView`. | Complete |
| 2 | Port scalar LaTeX generators. | Complete |
| 3 | Port LaTeX table generators. | Complete |
| 4 | Port figure data builders and plots. | Complete |
| 5 | Port PowerPoint slide generation for ISM molecules and PPD detections. | Pending |
| 6 | Add full-output orchestration and PDF-generation workflow. | Pending |

## Legacy Function Inventory

| Legacy function | Lines | Kind | Proposed destination | Migration status | Verification target |
| --- | ---: | --- | --- | --- | --- |
| `version` | 34-35 | Metadata helper | `astromol.__init__` or package metadata | Pending | Package version source |
| `updated` | 37-38 | Metadata helper | `astromol.census` or release metadata | Pending | Current database/release metadata |
| `make_all_plots` | 48-74 | Orchestrator | `astromol.manuscript` | Pending | Calls all modern plot products |
| `make_all_latex` | 77-112 | Orchestrator | `astromol.manuscript` | Pending | Calls all modern LaTeX products |
| `change_color` | 115-145 | Plot helper | `astromol.figures` | Pending | Unit-level visual/color check |
| `print_variables` | 147-237 | Inspection helper | `astromol.analysis` or CLI | Pending | Human-readable inventory output |
| `inspect` | 240-255 | Inspection helper | CLI/debug utility | Pending | Human-readable record dump |
| `summary` | 257-269 | Inspection helper | CLI/debug utility | Pending | Human-readable molecule summary |
| `cumu_det_plot` | 276-478 | Figure | `astromol.figures` | Complete | Verify 2021 cumulative series and preview figure |
| `cumu_det_natoms_plot` | 481-747 | Figure | `astromol.figures` | Complete | Verify 2021 atom-count series and preview figure |
| `det_per_year_per_atom` | 750-867 | Figure | `astromol.figures` | Complete | Verify 2021 bubble-rate figure |
| `facility_shares` | 870-1044 | Figure | `astromol.figures` | Complete | Verify 2021 facility-share pie grid |
| `scopes_by_year` | 1047-1202 | Figure | `astromol.figures` | Complete | Verify 2021 facility cumulative traces |
| `periodic_heatmap` | 1205-1436 | Figure | `astromol.figures` | Complete | Verify 2021 element counts and preview figure |
| `mass_by_wavelengths` | 1439-1673 | Figure | `astromol.figures` | Complete | Verify 2021 KDE and add 2026 box/points replacement |
| `mols_waves_by_atoms` | 1676-1818 | Figure | `astromol.figures` | Complete | Verify 2021 six-panel figure and add 2026 bubble-map replacement |
| `du_histogram` | 1821-1894 | Figure | `astromol.figures` | Complete | Verify 2021 histogram and add 2026 exact-value bar replacement |
| `type_pie_chart` | 1897-2114 | Figure | `astromol.figures` | Complete | Verify 2021 molecule-type ring chart |
| `source_pie_chart` | 2117-2311 | Figure | `astromol.figures` | Complete | Verify 2021 source-type ring chart |
| `indiv_source_pie_chart` | 2314-2500 | Figure | `astromol.figures` | Complete | Verify 2021 individual-source ring chart |
| `mol_type_by_source_type` | 2503-2697 | Figure | `astromol.figures` | Complete | Verify 2021 molecule-type-by-source pie grid |
| `du_by_source_type` | 2700-2850 | Figure | `astromol.figures` | Complete | Verify 2021 DU/source KDE figure |
| `rel_du_by_source_type` | 2853-2992 | Figure | `astromol.figures` | Complete | Verify 2021 relative-DU/source KDE figure |
| `mass_by_source_type` | 2995-3152 | Figure | `astromol.figures` | Complete | Verify 2021 mass/source KDE and add 2026 box/strip replacement |
| `waves_by_source_type` | 3155-3356 | Figure | `astromol.figures` | Complete | Verify 2021 wavelength/source pie grid |
| `kappas` | 3359-3417 | Figure | `astromol.figures` | Complete | Verify 2021 kappa histogram and 2026 production guide |
| `waves_pie_chart` | 3420-3556 | Figure | None | Skipped | Legacy helper not used in the 2021 paper |
| `make_ism_tables` | 3563-3720 | LaTeX table | `astromol.latex` | Complete | Reproduce 2021 ISM tables |
| `make_exgal_count` | 3722-3739 | LaTeX scalar | `astromol.latex` | Complete | Reproduce `nexgal.tex` value |
| `make_exgal_percent` | 3741-3758 | LaTeX scalar | `astromol.latex` | Complete | Reproduce `nexgalpercent.tex` value |
| `_make_exgal_sentence` | 3760-3783 | LaTeX prose fragment | `astromol.latex` | Deferred | Rebuild with table/prose generation |
| `make_det_count` | 3785-3802 | LaTeX scalar | `astromol.latex` | Complete | Reproduce `ndetects.tex` value |
| `make_elem_count` | 3805-3832 | LaTeX scalar | `astromol.latex` | Complete | Reproduce `nelems.tex` value |
| `make_ppd_count` | 3834-3859 | LaTeX scalar | `astromol.latex` | Complete | Reproduce `nppds.tex` value |
| `make_exo_count` | 3861-3885 | LaTeX scalar | `astromol.latex` | Complete | Reproduce `nexos.tex` value |
| `make_ices_count` | 3887-3911 | LaTeX scalar | `astromol.latex` | Complete | Reproduce `nices.tex` value |
| `make_ppd_isos_count` | 3913-3938 | LaTeX scalar | `astromol.latex` | Complete | Reproduce `nppdisos.tex` value |
| `make_exgal_table` | 3941-4140 | LaTeX table | `astromol.latex` | Complete | Reproduce 2021 exgal table |
| `make_ppd_table` | 4143-4301 | LaTeX table | `astromol.latex` | Complete | Reproduce 2021 PPD table |
| `make_exo_table` | 4304-4406 | LaTeX table | `astromol.latex` | Complete | Reproduce 2021 exoplanet table |
| `make_ice_table` | 4408-4523 | LaTeX table | `astromol.latex` | Complete | Reproduce 2021 ice table |
| `make_det_per_year_by_atoms_table` | 4525-4655 | LaTeX table | `astromol.latex` | Complete | Reproduce rate table |
| `make_facility_table` | 4657-4753 | LaTeX table | `astromol.latex` | Complete | Reproduce 2021 facility table |
| `make_source_table` | 4756-4860 | LaTeX table | `astromol.latex` | Complete | Reproduce 2021 source table |
| `make_rate_counts` | 4862-4913 | LaTeX scalar group | `astromol.latex` | Complete | Reproduce rate scalar files |
| `make_percent_radio` | 4915-4939 | LaTeX scalar | `astromol.latex` | Complete | Reproduce `radiopercent.tex` value |
| `make_scopes_count` | 4941-4964 | LaTeX scalar | `astromol.latex` | Complete | Reproduce `nscopes.tex` value |
| `make_percent_unsat` | 4966-4995 | LaTeX scalar | `astromol.latex` | Complete | Reproduce `unsatpercent.tex` value |
| `make_sat_list` | 4997-5026 | LaTeX scalar/list | `astromol.latex` | Complete | Reproduce `satlist.tex` content |
| `make_sat_count` | 5028-5054 | LaTeX scalar | `astromol.latex` | Complete | Reproduce `nsats.tex` value |
| `make_sat_percent` | 5056-5085 | LaTeX scalar | `astromol.latex` | Complete | Reproduce `satpercent.tex` value |
| `make_sfr_rad_percent` | 5087-5113 | LaTeX scalar | `astromol.latex` | Complete | Reproduce `sfr_rad_percent.tex` value |
| `make_dark_rad_percent` | 5115-5141 | LaTeX scalar | `astromol.latex` | Complete | Reproduce `dark_rad_percent.tex` value |
| `make_mols_slide` | 5147-5561 | PowerPoint slide | `astromol.slides` | Pending | Visual comparison to legacy slide |
| PPD detections slide | New | PowerPoint slide | `astromol.slides` | Pending | New product for PPD detection summary |

## Dependency Notes

Legacy imports include `numpy`, `matplotlib`, `periodictable`, `colour`,
`seaborn`, `scipy`, and `python-pptx`. The modern implementation should add
project dependency metadata before these become required package features.
Each dependency should be justified by the migrated code:

- Keep `numpy` for array/statistical helpers and for parity with legacy
  polynomial rate fits.
- Keep `matplotlib` for plotting unless a specific replacement is chosen.
- Keep `python-pptx` for the molecule slide product.
- Keep `scipy` for legacy KDE reproduction plots that still use
  `scipy.stats.gaussian_kde`.
- Re-evaluate `periodictable`, `colour`, and `seaborn` before adding package
  dependency metadata; avoid carrying dependencies that only supported legacy
  styling.

## Figure Migration Checkpoint

The database-backed figure phase has migrated every census-summary figure that
is included in the 2021 manuscript body:

- `cumulative_detections.pdf`
- `cumulative_by_atoms.pdf`
- `rate_by_atoms.pdf`
- `facility_shares.pdf`
- `scopes_by_year.pdf`
- `periodic_heatmap.pdf`
- `mass_by_wavelengths_kde.pdf`
- `mols_waves_by_atoms.pdf`
- `du_histogram.pdf`
- `type_pie_chart.pdf`
- `source_pie_chart.pdf`
- `mol_type_by_source_type.pdf`
- `du_by_source_type_kde.pdf`
- `relative_du_by_source_type_kde.pdf`
- `mass_by_source_type_kde.pdf`
- `waves_by_source_type.pdf`
- `kappas_histo.pdf`

Appendix/static explanatory figures from the 2021 arXiv bundle are not part of
the database figure migration. They can be reused as static assets or dropped
from the 2026 manuscript independently. The unused legacy `waves_pie_chart`
helper was also skipped because it was not included in the 2021 manuscript.

The next output-generation phase is slide generation, especially the legacy
ISM molecule slide and the new PPD detections slide.

## Completed View Verification

`test_census_view.py` verifies the shared view layer against the audited
historical table memberships:

- Census views exclude isotopologues by default and expose the isotope-expanded
  record with `include_isotopologues=True`.
- 2018 ISM/CSM: 204 secure molecules
- 2018 extragalactic: 63 secure, 65 with tentative rows
- 2018 exoplanet: 5 secure molecules
- 2018 PPD: 23 secure non-isotopologue molecules, 35 with isotopologues
- 2021 ISM/CSM: 240 secure molecules
- 2021 extragalactic: 73 secure, 75 with tentative rows
- 2021 exoplanet: 9 secure molecules
- 2021 PPD: 25 secure non-isotopologue molecules, 40 with isotopologues
- 2026/current ISM/CSM: 325 secure non-isotopologue molecules, 335 with
  isotopologues
- 2026/current exoplanet: 11 secure non-isotopologue molecules, 13 with
  isotopologues
- 2026/current PPD: 34 secure non-isotopologue molecules, 57 with isotopologues
- 2021 source/facility contribution counts
- 2026 census view and current view currently return identical secure
  accepted detection sets by context under both isotopologue settings

## Completed Scalar LaTeX Verification

`astromol.latex` now provides scalar fragment generators that consume a
`CensusView` and return filename-to-content mappings. `write_scalar_fragments`
can write those mappings to a manuscript input directory. The implementation
ports the count, percentage, saturated hydrocarbon, radical-source, telescope,
and detection-rate scalar outputs from the legacy module.

`test_latex_scalars.py` verifies the generated 2021 scalar fragments against
the audited refactor database view:

| Fragment | Generated value | Verified value | Result |
| --- | ---: | ---: | --- |
| `ndetects.tex` | 240 | 240 | Match |
| `nelems.tex` | 19 | 19 | Match |
| `nppds.tex` | 25 | 25 | Match |
| `nppdisos.tex` | 15 | 15 | Match |
| `nexgal.tex` | 73 | 73 | Match |
| `nexgalpercent.tex` | 30 | 30 | Match |
| `nexos.tex` | 9 | 9 | Match |
| `nices.tex` | 9 | 9 | Match |
| `radiopercent.tex` | 90 | 90 | Match |
| `nscopes.tex` | 46 | 46 | Match |
| `unsatpercent.tex` | 93 | 93 | Match |
| `satlist.tex` | `\ce{CH3Cl}, \ce{CH4}, \ce{CH3OH}, \ce{CH3SH}, \ce{CH3NH2}, \ce{CH3CH2OH}, \ce{CH3CH2SH}, \ce{CH3OCH3}, \ce{CH3OCH2OH}, and \ce{HOCH2CH2OH}` | Same | Match |
| `nsats.tex` | 10 | 10 | Match |
| `satpercent.tex` | 7 | 7 | Match |
| `sfr_rad_percent.tex` | 8 | 8 | Match |
| `dark_rad_percent.tex` | 25 | 25 | Match |
| `rate_since_1968.tex` | 3.9 | 3.9 | Match |
| `rate_since_2005.tex` | 6.0 | 6.0 | Match |

The legacy `_make_exgal_sentence` prose fragment is deferred to the table/prose
generation phase because its output depends on manuscript phrasing rather than
a standalone scalar value.

## Completed ISM Table Verification

`astromol.latex` now provides `ism_table_fragments` and `write_ism_tables`,
which port the legacy `make_ism_tables` output through `CensusView`. The
renderer supports `layout="legacy"` for audited 2021 reproduction and
`layout="balanced"` for production manuscript output. Table membership comes
from accepted non-isotopologue ISM/CSM detections in the view, while ordering
within each atom-count column uses the first accepted ISM/CSM detection date.
Molecule links use the compact `\molref{label}{formula}` macro instead of
repeating full `\hyperref...\ce...` markup in every table cell.

`test_latex_ism_tables.py` verifies the generated 2021 legacy ISM table
fragments:

| Check | Generated value | Verified value | Result |
| --- | ---: | ---: | --- |
| `ism_table_2-7.tex` molecule links | 186 | 186 | Match |
| `ism_table_8+.tex` molecule links | 54 | 54 | Match |
| Total linked 2021 ISM/CSM molecules | 240 | 240 | Match |
| `ism_table_2-7.tex` rendered rows | 23 | 23 | Match |
| `ism_table_8+.tex` rendered rows | 15 | 15 | Match |
| 2-7 atom column lengths | `21, 20, 23, 22, 16, 15, 16, 15, 23, 15` | Same | Match |
| 8+ atom column lengths | `15, 14, 6, 6, 5, 2, 3, 3` | Same | Match |
| Linked molecule-label set | 2021 accepted ISM/CSM labels | 2021 accepted ISM/CSM labels | Match |

The same test verifies the balanced 2026/current-oriented layout:

| Check | Generated value | Verified value | Result |
| --- | ---: | ---: | --- |
| Balanced table fragments | 3 | 3 | Match |
| Balanced table column counts | `7, 7, 6` | Non-increasing density | Match |
| Maximum columns per table | 7 | <= 7 | Match |
| Maximum rows per column | 23 | <= 23 | Match |
| Total linked 2026 ISM/CSM molecules | 325 | 325 | Match |
| Linked molecule-label set | 2026 accepted non-isotopologue ISM/CSM labels | 2026 accepted non-isotopologue ISM/CSM labels | Match |
| Duplicate linked labels | 0 | 0 | Match |

The main ISM molecule table intentionally excludes isotopologue records. The
balanced layout includes an explicit `13+ Atoms` category before PAHs and
fullerenes so large non-PAH/non-fullerene molecules remain represented as the
database grows. When a category is split across multiple columns, the column
header is rendered once with `\multicolumn` so the atom-count label remains
centered over the split columns, matching the 2021 manuscript style. A portrait
AASTeX preview of the current balanced tables was generated successfully during
the migration check; standalone preview warnings for unresolved molecule
hyperlinks are expected because the full manuscript section labels are not
present in the preview wrapper. The production layout packs earlier tables at
least as densely as later tables so the visual density stays constant or
decreases across table number; the current 2026/current-oriented output uses
seven, seven, and six columns across the three table fragments.

## Completed Extragalactic Table Verification

`astromol.latex` now provides `exgal_table_fragments` and `write_exgal_table`,
which port the legacy `make_exgal_table` output through `CensusView`. The table
includes secure extragalactic detections plus introduced tentative
extragalactic detections, and marks tentative rows with a dagger. It
intentionally excludes isotopologues, matching the standard census-table
policy. Atom-count groups are generated from the non-empty table bins so future
10-, 11-, or larger-atom exgal detections will be represented automatically.
Within each atom-count column, rows preserve database order to mirror the
legacy `all_molecules` filtering behavior. Observation references are numbered
by first appearance in the rendered table and then emitted as `\citet{...}`
notes below the table.

`test_latex_exgal_table.py` verifies the generated 2021 extragalactic table:

| Check | Generated value | Verified value | Result |
| --- | ---: | ---: | --- |
| Linked 2021 exgal molecules, including tentative | 75 | 75 | Match |
| Secure 2021 exgal molecules | 73 | 73 | Match |
| Tentative 2021 exgal rows | 2 | 2 | Match |
| Numbered citation entries | 52 | 52 | Match |
| First atom-group column lengths | `22, 18, 12, 9` | Same | Match |
| Second atom-group column lengths | `5, 5, 2, 1, 1` | Same | Match |
| Isotopologue rows | 0 | 0 | Match |
| Duplicate linked labels | 0 | 0 | Match |

The same test verifies the current 2026-oriented exgal table has 80 linked
non-isotopologue rows, two tentative rows, and 55 numbered citation entries.
A portrait AASTeX preview was generated successfully during the migration
check at `/private/tmp/astromol_exgal_table_preview/exgal_table_preview.pdf`.
Standalone preview warnings for unresolved molecule hyperlinks and citations
are expected because the full manuscript labels and bibliography are not
present in the preview wrapper; no overfull or underfull boxes were reported.

## Completed PPD Table Verification

`astromol.latex` now provides `ppd_table_fragments` and `write_ppd_table`,
which port the legacy `make_ppd_table` output through `CensusView`. Unlike the
main ISM/CSM and extragalactic molecule tables, the PPD table intentionally
includes isotopologue records. The table includes only secure PPD detections;
tentative and disputed detections remain excluded. Rows are ordered
parent-first, with detected isotopologues listed immediately after their parent
molecule where possible, matching the legacy nested `ppd_isos` behavior.
Atom-count groups are generated from non-empty bins so larger PPD detections
will be represented automatically as the inventory grows.

`test_latex_ppd_table.py` verifies the generated 2021 PPD table:

| Check | Generated value | Verified value | Result |
| --- | ---: | ---: | --- |
| Linked 2021 PPD molecules and isotopologues | 40 | 40 | Match |
| Non-isotopologue 2021 PPD molecules | 25 | 25 | Match |
| Isotopologue 2021 PPD rows | 15 | 15 | Match |
| Numbered citation entries | 38 | 38 | Match |
| Atom-group column lengths | `14, 16, 4, 4, 2` | Same | Match |
| Tentative/disputed rows | 0 | 0 | Match |
| Duplicate linked labels | 0 | 0 | Match |

The same test verifies the current 2026-oriented PPD table has 57 linked rows,
including 23 isotopologue rows, split into atom-count groups
`(2, 3, 4, 5, 6)` and `(7, 9, 12)`. A portrait AASTeX preview was generated
successfully during the migration check at
`/private/tmp/astromol_ppd_table_preview/ppd_table_preview.pdf`. Standalone
preview warnings for unresolved molecule hyperlinks and citations are expected
because the full manuscript labels and bibliography are not present in the
preview wrapper; no overfull or underfull boxes were reported.

## Completed Exoplanet Table Verification

`astromol.latex` now provides `exoplanet_table_fragments` and
`write_exoplanet_table`, which port the legacy `make_exo_table` output through
`CensusView`. The standard exoplanet table is isotope-free, matching the
default census-table policy and the legacy table behavior. An explicit
`include_isotopologues=True` option is available for expanded exoplanet views
when needed. The table includes only secure exoplanet-atmosphere detections;
tentative and disputed detections remain excluded. Rows preserve database order
to mirror the legacy `all_molecules` filtering behavior. Observation references
are numbered by first appearance in the rendered table and emitted as
`\citet{...}` notes below the table.

`test_latex_exoplanet_table.py` verifies the generated 2021 exoplanet table:

| Check | Generated value | Verified value | Result |
| --- | ---: | ---: | --- |
| Linked 2021 exoplanet molecules | 9 | 9 | Match |
| Isotopologue rows in standard table | 0 | 0 | Match |
| Numbered citation entries | 19 | 19 | Match |
| Tentative/disputed rows | 0 | 0 | Match |
| Duplicate linked labels | 0 | 0 | Match |

The same test verifies the current 2026-oriented standard exoplanet table has
11 non-isotopologue rows and that the optional isotope-expanded output has 13
rows, including `^{13}CO` and `CH3D`. A portrait AASTeX preview was generated
successfully during the migration check at
`/private/tmp/astromol_exoplanet_table_preview/exo_table_preview.pdf`.
Standalone preview warnings for unresolved molecule hyperlinks and citations
are expected because the full manuscript labels and bibliography are not
present in the preview wrapper; no overfull or underfull boxes were reported.

## Completed Ice Table Verification

`astromol.latex` now provides `ice_table_fragments` and `write_ice_table`,
which port the legacy `make_ice_table` output through `CensusView`. The
standard ice table includes secure and tentative detections, excludes
isotopologue records, and marks tentative rows with a dagger. Pass
`include_tentative=False` for an explicit secure-only output.

One legacy nuance is intentionally handled by status rather than a manual table
insertion. The 2021 script inserted `OCN-` directly into the ice table, but the
modern database records the OCN- ice claim as tentative. It is therefore
included in standard ice tables but marked with a dagger.

`test_latex_ice_table.py` verifies the generated 2021 ice table:

| Check | Generated value | Verified value | Result |
| --- | ---: | ---: | --- |
| Linked 2021 ice molecules | 10 | 10 | Match |
| Isotopologue rows in standard table | 0 | 0 | Match |
| Numbered citation entries | 11 | 11 | Match |
| Tentative rows in standard table | 1 | 1 | Match |
| Duplicate linked labels | 0 | 0 | Match |
| OCN- in standard table | 1 | 1 | Match |

The same test verifies the current 2026-oriented standard ice table has 15
non-isotopologue rows, including `OCN-` with a dagger marker, and 12 numbered
citation entries. The secure-only output has 14 rows and excludes `OCN-`.
A portrait AASTeX preview was generated successfully during the migration check
at `/private/tmp/astromol_ice_table_preview/ice_table_preview.pdf`.
Standalone preview warnings for unresolved molecule hyperlinks and citations
are expected because the full manuscript labels and bibliography are not
present in the preview wrapper; no overfull or underfull boxes were reported.

## Completed Detection-Rate-By-Atoms Table Verification

`astromol.latex` now provides `rate_by_atoms_table_fragments` and
`write_rate_by_atoms_table`, which port the legacy
`make_det_per_year_by_atoms_table` output through `CensusView`. The table uses
secure, non-isotopologue ISM/CSM first-detection years. Two legacy assumptions
were modernized deliberately: the large normal-molecule bin is now a real `13+`
category that excludes PAHs and fullerenes, and the caption no longer refers to
`scipy.stats.linregress`; fits are calculated with NumPy. The legacy table
header said `R^2`, but the generated values were SciPy's Pearson `r_value`.
The modern table instead reports true `R^2` values. The 2026 manuscript should
note that the previous paper's rate table mistakenly displayed Pearson `R`
values under an `R^2` heading.

`test_latex_rate_by_atoms_table.py` verifies the generated 2021 table:

| Check | Generated value | Verified value | Result |
| --- | ---: | ---: | --- |
| Fitted rows | 13 | 13 | Match |
| 2-atom rate and R2 | `0.68`, `0.98` | Same | Match |
| 12-atom rate and R2 | `0.16`, `0.84` | Same | Match |
| 13+ rate and R2 | `0.30`, `0.60` | Same | Match |
| Fullerene rate and R2 | `0.09`, `0.57` | Same | Match |
| PAH row in 2021 table | omitted | omitted | Match |

The 2021 PAH row is omitted because all 2021 PAH detections occur in the onset
year, leaving only one point for that census-bounded fit. The current
2026-oriented output has both PAH and fullerene rows; the `13+` row includes
normal non-PAH/non-fullerene molecules with 13 or more atoms. Portrait AASTeX
previews were generated successfully during the migration check at
`/private/tmp/astromol_rate_by_atoms_2021_preview/rates_by_atoms_preview.pdf`
and
`/private/tmp/astromol_rate_by_atoms_2026_preview/rates_by_atoms_preview.pdf`.
Standalone preview warnings for unresolved figure references are expected
because the full manuscript labels are not present in the preview wrappers.

## Completed Facility Table Verification

`astromol.latex` now provides `facility_table_fragments` and
`write_facility_table`, which port the legacy `make_facility_table` output
through `CensusView`. The table uses secure, non-isotopologue ISM/CSM
detections and counts each observing facility listed on those detection
records. Facilities with no detections in the selected view are omitted. Rows
are sorted by descending count with alphabetical tie-breaking and rendered as
two side-by-side facility/count column pairs.

`test_latex_facility_table.py` verifies the generated 2021 facility table:

| Check | Generated value | Verified value | Result |
| --- | ---: | ---: | --- |
| Facilities with detections | 46 | 46 | Match |
| Total facility contributions | 304 | 304 | Match |
| IRAM 30-m count | 64 | 64 | Match |
| NRAO 36-ft count | 33 | 33 | Match |
| GBT 100-m count | 28 | 28 | Match |
| NRAO/ARO 12-m count | 27 | 27 | Match |
| Yebes 40-m count | 19 | 19 | Match |

The current 2026-oriented output has 47 facilities with detections and 415
facility contributions. Portrait AASTeX previews were generated successfully
during the migration check at
`/private/tmp/astromol_facility_table_2021_preview/facilities_preview.pdf` and
`/private/tmp/astromol_facility_table_2026_preview/facilities_preview.pdf`.

## Completed Source Table Verification

`astromol.latex` now provides `source_table_fragments` and
`write_source_table`, which port the legacy `make_source_table` output through
`CensusView`. The table uses secure, non-isotopologue ISM/CSM detections and
counts each source listed on those detection records. Diffuse-cloud
line-of-sight sources are consolidated into a single `Diffuse Cloud` row, and
source-region grouping such as Sgr B2 is inherited from the curated source
records. Sources with no detections in the selected view are omitted. Rows are
sorted by descending count with alphabetical tie-breaking and rendered as two
side-by-side source/count column pairs.

`test_latex_source_table.py` verifies the generated 2021 source table:

| Check | Generated value | Verified value | Result |
| --- | ---: | ---: | --- |
| Sources with detections | 38 | 38 | Match |
| Total source contributions | 332 | 332 | Match |
| Sgr B2 count | 68 | 68 | Match |
| TMC-1 count | 57 | 57 | Match |
| IRC+10216 count | 55 | 55 | Match |
| Diffuse Cloud count | 42 | 42 | Match |
| Orion count | 24 | 24 | Match |

The current 2026-oriented output has 42 sources with detections and 424 source
contributions. Portrait AASTeX previews were generated successfully during the
migration check at
`/private/tmp/astromol_source_table_2021_preview/source_preview.pdf` and
`/private/tmp/astromol_source_table_2026_preview/source_preview.pdf`.

## Completed Cumulative-Detections Figure Verification

`astromol.figures` now provides `cumulative_detection_data`,
`plot_cumulative_detections`, and `write_cumulative_detections_plot`, which
port the legacy `cumu_det_plot` output through `CensusView`. The data layer is
separate from matplotlib rendering so the cumulative series and fitted rates
can be tested directly. By default, the figure uses secure, non-isotopologue
ISM/CSM molecule detections. The 2021 verification view ends at the 2021 census
boundary; the current 2026-oriented preview ends at the 2026 census boundary.

The linear trend ranges are census-aware and use closed date ranges:
2018 and 2021 legacy views use `1968-2005` plus `2005-<census>`, the 2026
census view uses `1968-2005`, `2005-2021`, and `2021-2026`, and the live
current view uses `1968-2005`, `2005-2021`, and `2021-Present`. Facility
commissioning annotations are retained for the 2018 and 2021 legacy views but
are omitted from 2026/current outputs. Trend annotations and total-count labels
are anchored in axes-fraction coordinates and the y-axis has explicit headroom,
so the cumulative curve does not crowd the figure text as the database grows.

`test_figures_cumulative_detections.py` verifies the generated 2021 cumulative
detection data:

| Check | Generated value | Verified value | Result |
| --- | ---: | ---: | --- |
| Start year | 1937 | 1937 | Match |
| End year | 2021 | 2021 | Match |
| Final cumulative count | 240 | 240 | Match |
| 1968 count | 5 | 5 | Match |
| 2004 count | 130 | 130 | Match |
| 2019 count | 213 | 213 | Match |
| 2020 count | 219 | 219 | Match |
| 1968-2005 rate | 3.5 | 3.5 | Match |
| 2005-2021 rate | 6.0 | 6.0 | Match |

The current 2026-oriented output has a final cumulative count of 325. Its
annotated rates are 3.5 detections/year for 1968-2005, 6.4 detections/year for
2005-2021, and 14.2 detections/year for 2021-2026. Preview PDFs were generated
successfully during the migration check at
`/private/tmp/astromol_cumulative_detections_2021_preview/cumulative_detections.pdf`
and
`/private/tmp/astromol_cumulative_detections_2026_preview/cumulative_detections.pdf`.

Initial assessment: the legacy figure remains useful as a high-level narrative
plot, but facility annotations are too crowded for the 2026/current versions.
A complementary figure or inset focused on the post-2018 acceleration may be
more informative than adding more annotations to the cumulative plot itself.

## Completed Cumulative-By-Atoms Figure Verification

`astromol.figures` now provides `cumulative_by_atoms_data`,
`plot_cumulative_by_atoms`, and `write_cumulative_by_atoms_plot`, porting the
legacy `cumu_det_natoms_plot` output through `CensusView`. The data layer uses
the same secure, non-isotopologue ISM/CSM first-detection records as the main
cumulative plot and groups them into 2-12 atom traces, a true non-PAH/non-
fullerene `13+` trace, and separate PAH and fullerene traces.

The same data layer also powers a promoted companion product:
`plot_stacked_cumulative_by_atoms` and `write_stacked_cumulative_by_atoms_plot`.
This stacked view is not a 2021 legacy reproduction target; it is a
composition-focused manuscript candidate for showing how atom-count classes
contribute to total inventory growth. It uses the same colors and the same
in-plot legend placement as the legacy atom-count line plot. The line and
stacked outputs use the same figure size, fixed axes rectangle, x-limits, and
legend coordinates so they can be placed side-by-side as matched panels in a
single LaTeX figure. Common figure size, tick, axes, and text styling now live
in shared helpers in `astromol.figures` so future plots can keep a coherent
manuscript style.

Atom-count figure colors are view-aware. Legacy 2018/2021 reproduction views
keep the original palette for visual parity, while 2026/current and future
views use a solid color-blind-friendly palette. Dashed lines and hatches were
tested as accessibility aids but rejected because they made these already dense
figures harder to read.

For the 2021 verification view, the modern `13+` definition does not change
the legacy count because there are no non-PAH/non-fullerene molecules with more
than 13 atoms. For the 2026 view it does matter: `mol:c-C6H5CCH` has 14 atoms
and is included in the modern `13+` trace.

`test_figures_cumulative_by_atoms.py` verifies the generated 2021 final
category counts:

| Category | Count |
| --- | ---: |
| 2 atoms | 41 |
| 3 atoms | 45 |
| 4 atoms | 31 |
| 5 atoms | 31 |
| 6 atoms | 23 |
| 7 atoms | 15 |
| 8 atoms | 15 |
| 9 atoms | 14 |
| 10 atoms | 6 |
| 11 atoms | 6 |
| 12 atoms | 5 |
| 13+ atoms | 2 |
| Fullerenes | 3 |
| PAHs | 3 |

Preview PDFs were generated successfully at
`/private/tmp/astromol_cumulative_by_atoms_2021_preview/cumulative_by_atoms.pdf`
and
`/private/tmp/astromol_cumulative_by_atoms_2026_preview/cumulative_by_atoms.pdf`.
Stacked companion previews were generated at
`/private/tmp/astromol_stacked_cumulative_by_atoms_2021_preview/stacked_cumulative_by_atoms.pdf`
and
`/private/tmp/astromol_stacked_cumulative_by_atoms_2026_preview/stacked_cumulative_by_atoms.pdf`.

The atom-count data layer also now supports a promoted rolling-rate heatmap:
`rolling_rate_by_atoms_heatmap_data`,
`plot_rolling_rate_by_atoms_heatmap`, and
`write_rolling_rate_by_atoms_heatmap`. This figure is not a legacy
reproduction target. It was selected from exploratory rolling-window and
time-binned alternatives because it most directly shows when each atom-count
class accelerated or quieted while remaining readable in a single-column
layout. The production version uses a 10-year trailing rate, the `cividis`
sequential palette, fixed 0-2 detections/year color scaling, no row separator
lines, and short row labels. The fullerene row is labeled `Fuller` to identify
the class without implying the row is only C60.

`test_figures_cumulative_by_atoms.py` verifies the heatmap data shape, labels,
2026 final-window rates, fixed canvas geometry, absence of row separator
lines, and file writing. Preview PDFs were generated successfully at
`/private/tmp/astromol_rolling_rate_heatmap_2021_preview/rolling_rate_heatmap.pdf`,
`/private/tmp/astromol_rolling_rate_heatmap_2026_preview/rolling_rate_heatmap.pdf`,
and
`/private/tmp/astromol_rolling_rate_heatmap_current_preview/rolling_rate_heatmap.pdf`.

## Completed Detection-Rate-By-Atoms Figure Verification

`astromol.figures` now provides `detection_rate_by_atoms_data`,
`plot_detection_rate_by_atoms`, and `write_detection_rate_by_atoms_plot`,
porting the legacy `det_per_year_per_atom` bubble plot through `CensusView`.
The figure uses secure, non-isotopologue ISM/CSM first-detection records. For
each atom-count category, it divides the number of detections by the number of
years from that category's first detection through the view end year. Marker
area remains proportional to the total number of detections in the category.
The legacy single-view writer remains available for reproducing the 2021
artifact.

The promoted production replacement is
`plot_detection_rate_by_atoms_comparison` and
`write_detection_rate_by_atoms_comparison_plot`. It overlays the target
census/current view on a 2021 baseline: the target view is rendered as filled
blue circles, and the 2021 values are rendered behind it as unfilled hatched
black circles. No horizontal offset is used. Unlike the legacy implementation,
the comparison plot uses one absolute marker-area scale for both datasets
(`area = count * constant`), so marker area can be compared directly between
2021 and the target view. This scaling change must be described in the 2026
manuscript if the figure is compared to the 2021 version.

The 2021 reproduction intentionally preserves a legacy plotting convention:
13+, PAH, and fullerene categories are included in the data and marker-scaling
calculation, but the displayed x-axis shows only 2-12 atom categories. The
modern single-view implementation draws only visible points so hidden
categories cannot bleed into the right edge of 2026/current previews as marker
sizes change. The production comparison plot also keeps the 2-12 atom x-axis,
because the unclipped 2026 view makes PAHs visually dominate the figure and
obscures the normal atom-count progression.

`test_figures_detection_rate_by_atoms.py` verifies the generated 2021 category
counts, first-detection years, rates, x-axis range, marker rendering, and file
writing. Key 2021 visible rates are:

| Category | Count | First year | Detections/year |
| --- | ---: | ---: | ---: |
| 2 atoms | 41 | 1937 | 0.48 |
| 3 atoms | 45 | 1969 | 0.85 |
| 4 atoms | 31 | 1968 | 0.57 |
| 5 atoms | 31 | 1971 | 0.61 |
| 6 atoms | 23 | 1970 | 0.44 |
| 7 atoms | 15 | 1973 | 0.31 |
| 8 atoms | 15 | 1975 | 0.32 |
| 9 atoms | 14 | 1974 | 0.29 |
| 10 atoms | 6 | 1987 | 0.17 |
| 11 atoms | 6 | 1978 | 0.14 |
| 12 atoms | 5 | 2001 | 0.24 |

Legacy single-view preview PDFs were generated successfully at
`/private/tmp/astromol_detection_rate_by_atoms_2021_preview/rate_by_atoms.pdf`,
`/private/tmp/astromol_detection_rate_by_atoms_2026_preview/rate_by_atoms.pdf`,
and
`/private/tmp/astromol_detection_rate_by_atoms_current_preview/rate_by_atoms.pdf`.
Production comparison preview PDFs were generated successfully at
`/private/tmp/astromol_detection_rate_by_atoms_comparison_2026_preview/rate_by_atoms_comparison.pdf`
and
`/private/tmp/astromol_detection_rate_by_atoms_comparison_current_preview/rate_by_atoms_comparison.pdf`.
The figure remains visually simple and does not present a color-accessibility
problem. The hatched 2021 background markers add non-color redundancy, while
the current/2026 values retain the package blue.

## Completed Facility-Shares Figure Verification

`astromol.figures` now provides `facility_share_data`,
`plot_facility_shares`, `write_facility_shares_plot`,
`plot_facility_share_bars`, and `write_facility_share_bars_plot`. The pie-grid
function ports the legacy `facility_shares` figure through `CensusView`; the
bar chart is the preferred production-facing view because it makes rank,
relative magnitude, and denominators easier to compare. The data layer uses
secure, non-isotopologue ISM/CSM first-detection records. A facility's share is
computed as its first-detection contribution count divided by all first
detections that occurred during that facility's operational window. Active
facilities are shown in the project default blue, while facilities
decommissioned by the selected view boundary are shown in the legacy pastel
red (`#F87070`). The legacy pie chart intentionally keeps its square 3x3
aspect ratio instead of the standard single-panel figure aspect ratio.

The 2021 reproduction preserves the published facility set:

| Facility | Window | Contributions | Window detections | Share |
| --- | --- | ---: | ---: | ---: |
| NRAO 36-ft | 1967-1984 | 33 | 56 | 58% |
| IRAM 30-m | 1984-2021 | 64 | 184 | 34% |
| GBT 100-m | 2004-2021 | 28 | 115 | 24% |
| Herschel | 2009-2013 | 7 | 33 | 21% |
| Yebes 40-m | 2007-2021 | 19 | 102 | 18% |
| NRAO/ARO 12-m | 1984-2021 | 27 | 184 | 14% |
| Bell 7-m | 1976-1992 | 8 | 59 | 13% |
| ALMA | 2011-2021 | 7 | 76 | 9% |
| NRAO 140-ft | 1965-2008 | 13 | 146 | 8% |

The modern curated 2021 top-nine selection would include Nobeyama 45-m and
drop Herschel, because current first-detection links assign 15 first
detections to Nobeyama. This is recorded as a data/legacy-artifact difference;
the default 2021 view preserves the published figure, while 2026/current views
use the curated top-nine selection.

`test_figures_facility_shares.py` verifies the 2021 facility set, percentage
values, active/inactive classification, the Nobeyama/Herschel selection
difference, pie-grid rendering, bar-chart rendering, and file writing. Preview
PDFs for the accepted pie-grid reproduction were generated successfully at
`/private/tmp/astromol_facility_shares_2021_preview/facility_shares.pdf`,
`/private/tmp/astromol_facility_shares_2026_preview/facility_shares.pdf`, and
`/private/tmp/astromol_facility_shares_current_preview/facility_shares.pdf`.
Production bar-chart preview PDFs were generated successfully at
`/private/tmp/astromol_facility_share_bars_production_preview/facility_share_bars_2021.pdf`,
`/private/tmp/astromol_facility_share_bars_production_preview/facility_share_bars_2026.pdf`,
and
`/private/tmp/astromol_facility_share_bars_production_preview/facility_share_bars_current.pdf`.

## Completed Scopes-By-Year Figure Verification

`astromol.figures` now provides `scopes_by_year_data`,
`plot_scopes_by_year`, and `write_scopes_by_year_plot`, porting the legacy
`scopes_by_year` cumulative facility-contribution figure through
`CensusView`. The data layer uses secure, non-isotopologue ISM/CSM first
detections and counts a contribution for each telescope listed on a first
detection. Facilities are included when their contribution count is at least
the configurable threshold, matching the legacy default of 10.

The 2021 reproduction matches the legacy facility set, cumulative final
counts, and rounded linear-rate annotations:

| Facility | Final count | Rate | Fit window |
| --- | ---: | ---: | --- |
| NRAO 36-ft | 33 | 2.2/yr | 1967-1985 |
| IRAM 30-m | 64 | 1.5/yr | 1984-2021 |
| GBT 100-m | 28 | 1.1/yr | 2004-2021 |
| Nobeyama 45-m | 15 | 1.1/yr | 1982-1997 |
| NRAO/ARO 12-m | 27 | 0.8/yr | 1984-2021 |
| Yebes 40-m | 19 | 0.7/yr | 2007-2021 |
| NRAO 140-ft | 13 | 0.5/yr | 1965-1993 |

The 2026 preview crosses the same threshold for ALMA and shows Yebes 40-m as
the highest-rate facility in the current database state. The accepted modern
production style keeps legacy reproductions available while using the 2026
figure palette, emphasizing Yebes 40-m in MIT Red, and drawing 50%-opacity
dotted post-last-contribution tails for NRAO 36-ft, Nobeyama 45-m, and NRAO
140-ft. The dotted tails indicate no later credited first detections in the
selected view, not necessarily telescope inactivity.

A complementary facility-era stacked bar chart was previewed and rejected; it
is not retained as a production output.
`test_figures_scopes_by_year.py` verifies the 2021 final counts, rounded
rates, fit windows, modern Yebes highlighting, dormant dotted tails, basic
plot styling, and file writing.

## Completed Periodic-Heatmap Figure Verification

`astromol.figures` now provides `periodic_heatmap_data`,
`plot_periodic_heatmap`, and `write_periodic_heatmap`, porting the legacy
`periodic_heatmap` figure through `CensusView`. The data layer uses secure,
non-isotopologue ISM/CSM molecules by default.

The 2021 data reproduce the arXiv scalar files: 240 ISM/CSM molecules and 19
detected elements. The largest element counts are C=188, H=173, N=92, O=75,
S=30, and Si=13. Preview PDFs were generated successfully at
`/private/tmp/astromol_periodic_heatmap_preview/periodic_heatmap_2021.pdf` and
`/private/tmp/astromol_periodic_heatmap_preview/periodic_heatmap_2026.pdf`.
`test_figures_periodic_heatmap.py` verifies the scalar counts, major element
counts, periodic-cell placement, basic plot construction, and file writing.
The faithful yellow-to-red legacy heatmap was accepted as the production
version. `cividis` and `viridis` palette variants were previewed and rejected.

## Completed Mass-By-Wavelength Figure Migration

`astromol.figures` now provides `mass_by_wavelength_data`,
`plot_mass_by_wavelength`, and `write_mass_by_wavelength_plot`, porting the
legacy `mass_by_wavelengths_kde` figure through `CensusView`. The data layer
uses secure, non-isotopologue ISM/CSM first detections by default and the
legacy KDE bandwidth convention (`bw = 0.5`).

The 2021 generated counts match the published figure labels for cm, mm,
sub-mm, and IR, but not for UV/Vis:

| Wavelength group | Published 2021 label | Current generated 2021 count |
| --- | ---: | ---: |
| cm | 89 | 89 |
| mm | 137 | 137 |
| sub-mm | 13 | 13 |
| IR | 18 | 18 |
| UV/Vis | 5 | 6 |

The IR count confirms that the published figure included the three fullerene
IR detections despite the caption text saying the figure excepted fullerenes.
The UV/Vis difference is a legacy plotting artifact: the old function merged
the UV and visible lists by molecular mass rather than by molecule identity,
and the legacy mass calculation ignored charge, causing CH+ to be dropped
after CH.
Preview PDFs were generated successfully at
`/private/tmp/astromol_mass_by_wavelength_preview/mass_by_wavelengths_kde_2021.pdf`
and
`/private/tmp/astromol_mass_by_wavelength_preview/mass_by_wavelengths_kde_2026.pdf`.
`test_figures_mass_by_wavelength.py` verifies the current database-derived
counts, basic plot construction, and file writing while documenting the known
2021 count mismatch.

The KDE path is retained for historical comparison with the 2021 census. For
the 2026 manuscript, the production replacement is
`plot_mass_by_wavelength_boxplot`/`write_mass_by_wavelength_boxplot`, a
horizontal box/whisker plus jittered-detection view. That view uses
`include_fullerenes=False` to match the intended 2021 caption language and to
avoid letting the fullerene IR detections dominate a figure meant to compare
mass distributions by discovery wavelength. Boxes show the interquartile
range, whiskers span the 10th-90th percentiles, individual points show the
underlying detections, and the 80 amu marker tracks the threshold discussed in
the manuscript text.

## Completed Molecules-By-Wavelength/Atom-Count Figure Migration

`astromol.figures` now provides `molecules_by_wavelength_atoms_data`,
`plot_molecules_by_wavelength_atoms`, and
`write_molecules_by_wavelength_atoms_plot`, porting the legacy six-panel
`mols_waves_by_atoms` KDE/histogram figure through `CensusView`. The data
layer uses secure, non-isotopologue ISM/CSM first detections by default and
excludes fullerenes by default for this figure.

The 2021 verification data reproduce the expected molecule count of 237 and
the wavelength-category counts cm=89, mm=137, sub-mm=13, IR=15, Vis=2, and
UV=6. During verification, wavelength mismatches against the legacy record
were traced to curated database corrections rather than plotting logic.
`test_figures_molecules_by_wavelength_atoms.py` verifies those counts, selected
matrix cells, legacy plot construction, and output writing.

The legacy six-panel generator is retained for historical reproduction. For
the 2026 manuscript, the accepted production replacement is
`plot_molecules_by_wavelength_atoms_bubble_heatmap`/
`write_molecules_by_wavelength_atoms_bubble_heatmap`, a count bubble map using
the same yellow-to-red count palette as the periodic-table heatmap. The
accepted view shows absolute counts of secure non-isotopologue non-fullerene
ISM/CSM first detections by wavelength category and atom-count bin. A
row-percentage preview mode was implemented for comparison but was not selected
as the manuscript figure because it obscured the absolute population
differences that the figure is intended to communicate.

## Completed Degree-Of-Unsaturation Figure Migration

`astromol.figures` now provides `du_histogram_data`,
`plot_du_histogram`, and `write_du_histogram`, porting the legacy
`du_histogram` through `CensusView`. The data layer uses secure,
non-isotopologue ISM/CSM molecules by default, excludes fullerenes by default,
and includes molecules whose formulas contain only H, D, N, C, Cl, F, S, and
O. This preserves the legacy DU-compatible formula domain.

The 2021 verification data contain 196 DU-compatible molecules. The legacy
histogram reproduction matches the 2021 visual artifact, including the
half-DU-width bins and the HC11N high-DU annotation. Four small protonated
hydride/non-carbon species (`H2Cl+`, `H3+`, `H3O+`, and `NH3D+`) evaluate to
DU = -0.5 under the formula expression; this is a domain artifact of applying
a neutral/organic DU expression outside its chemically meaningful range, not a
physical negative unsaturation.

For the 2026 manuscript, the accepted production replacement is
`plot_du_bar_chart`/`write_du_bar_chart`, an exact-value DU bar chart rather
than a histogram. The bar chart excludes negative-DU domain artifacts by
default, restores the saturated-species label at DU = 0, and labels HC11N and
the DU = 14 cyanopyrene isomers (`1-/2-/4-C16H9CN`) without per-bar numeric
labels. Several more radical alternatives, including cumulative exceedance,
DU-vs-atom-count scatter, categorical composition, and ranked DU strip views,
were previewed and rejected.

`test_figures_du_histogram.py` verifies the 2021 DU counts, saturated and
unsaturated counts, the negative-DU domain artifact, the DU = 14 cyanopyrene
labels, legacy histogram construction, exact-value bar chart construction, and
file writing.

## Completed Molecule-Type Figure Migration

`astromol.figures` now provides `molecule_type_data`,
`plot_type_pie_chart`, and `write_type_pie_chart`, porting the legacy
`type_pie_chart` concentric molecule-type ring figure through `CensusView`.
The data layer uses secure, non-isotopologue ISM/CSM molecules by default.
Categories are counted independently, so a molecule can contribute to more
than one category; this preserves the manuscript convention for neutral
radicals, cyclic ions, PAHs, and fullerenes.

The 2021 verification data contain 240 molecules and reproduce the legacy
category counts: 204 neutral, 54 radical, 30 cation, 19 cyclic, 6 anion, 3
fullerene, and 3 PAH. The 2026/current preview contains 325 molecules with
270 neutral, 71 radical, 47 cation, 31 cyclic, 8 anion, 3 fullerene, and 9 PAH.
The faithful legacy-style ring chart was accepted without a separate modern
replacement. Production rendering sorts rings by descending category count
while retaining stable category colors; `order="legacy"` remains available for
historical visual comparison.

`test_figures_type_pie_chart.py` verifies the 2021 and 2026 category counts,
rounded percentages, production ring ordering, plot labels, and file writing.

## Completed Source-Type Figure Migration

`astromol.figures` now provides `source_type_data`,
`plot_source_pie_chart`, and `write_source_pie_chart`, porting the legacy
`source_pie_chart` concentric source-type ring figure through `CensusView`.
The data layer uses secure, non-isotopologue ISM/CSM first detections by
default. Each molecule can receive credit for multiple generalized source
categories if its first detection lists sources in multiple categories, but it
is credited at most once per category.

The 2021 verification data contain 240 first detections and reproduce the
legacy category counts and percentages: SFR 87 (36.2%), dark cloud 67 (27.9%),
carbon star 58 (24.2%), other 32 (13.3%), and LOS/diffuse cloud 24 (10.0%).
The 2026/current preview contains 325 first detections with dark cloud 117
(36.0%), SFR 92 (28.3%), carbon star 71 (21.8%), other 50 (15.4%), and diffuse
cloud 24 (7.4%). Production rendering sorts rings by descending category count
while retaining stable category colors; this moves dark clouds to the outer
ring in the current 2026 view. `order="legacy"` remains available for
historical visual comparison.

`test_figures_source_pie_chart.py` verifies the 2021 and 2026 category counts,
rounded percentages, production ring ordering, plot labels, and file writing.

## Completed Individual-Source Figure Migration

`astromol.figures` now provides `individual_source_data`,
`plot_individual_source_pie_chart`, and
`write_individual_source_pie_chart`, porting the legacy
`indiv_source_pie_chart` concentric source ring figure through `CensusView`.
The data layer uses secure, non-isotopologue ISM/CSM first detections by
default.

This figure preserves the legacy source-contribution convention rather than
counting unique molecules per source. Named first-detection source
contributions are credited to Sgr B2, TMC-1, IRC+10216, and Orion; every other
listed first-detection source contribution is credited to Other. The
denominator remains the number of first-detected molecules, so percentages can
sum to more than 100% when multiple sources are listed for one first detection.

The 2021 verification data contain 240 first detections and reproduce the
legacy counts and percentages: Other 128 (53.3%), Sgr B2 68 (28.3%), TMC-1 57
(23.8%), IRC+10216 55 (22.9%), and Orion 24 (10.0%). The 2026/current preview
contains 325 first detections with Other 155 (47.7%), TMC-1 106 (32.6%), Sgr B2
70 (21.5%), IRC+10216 68 (20.9%), and Orion 25 (7.7%). Production rendering
sorts rings by descending count while retaining stable category colors;
`order="legacy"` remains available for historical visual comparison.

`test_figures_individual_source_pie_chart.py` verifies the 2021 and 2026
source contribution counts, rounded percentages, production ring ordering, plot
labels, and file writing.

## Completed Molecule-Type-By-Source Figure Migration

`astromol.figures` now provides `molecule_type_by_source_type_data`,
`plot_molecule_type_by_source_type`, and
`write_molecule_type_by_source_type`, porting the legacy
`mol_type_by_source_type` four-panel pie-grid figure through `CensusView`. The
data layer uses secure, non-isotopologue ISM/CSM first detections by default.

The source categories are the generalized first-detection source categories
used by the source-type figure: carbon star, dark cloud, diffuse cloud, and
SFR. A molecule is credited at most once per source category. Within each
source category, molecule-type classes are counted independently, so neutral
radicals, cyclic ions, and other overlapping classes can contribute to more
than one wedge.

The 2021 verification data reproduce the legacy wedge counts:
carbon star 6 anion, 3 cyclic, 52 neutral, 19 radical; dark cloud 2 anion, 9
cation, 10 cyclic, 56 neutral, 17 radical; LOS/diffuse cloud 6 cation, 18
neutral, 9 radical; and SFR 13 cation, 4 cyclic, 74 neutral, 7 radical. The
2026/current preview updates those counts while retaining the same visual
structure and stable molecule-type colors.

For the 2026 manuscript, the accepted production replacement is
`plot_molecule_type_by_source_enrichment_matrix`/
`write_molecule_type_by_source_enrichment_matrix`. This matrix keeps the same
source/type question but shows fractional enrichment relative to the overall
secure ISM/CSM molecule-type mix. The cell label gives the enrichment factor
and raw count; the color scale is log-centered on `1x`, with gray for
depletion and astromol blue for enrichment. This was chosen because it conveys
the same scientific point as the legacy pie grid more directly: which source
types are enriched or depleted in each molecule class.

`test_figures_molecule_type_by_source_type.py` verifies the 2021 and 2026
source-by-type counts, source denominators, overall type baseline counts,
legacy LOS label override, pie-grid labels, enrichment-matrix labels, legend,
and file writing.

## Completed DU-By-Source-Type Figure Migration

`astromol.figures` now provides `du_by_source_type_data`,
`plot_du_by_source_type`, and `write_du_by_source_type`, porting the legacy
`du_by_source_type` KDE figure through `CensusView`. It also provides
`plot_du_by_source_type_boxplot` and `write_du_by_source_type_boxplot` as the
preferred 2026 production replacement. The data layer uses secure,
non-isotopologue ISM/CSM first detections by default. Each molecule is credited
at most once per generalized first-detection source category, and fullerenes
are excluded by default so their extreme DU values do not dominate the
comparison.

The 2021 verification data reproduce the legacy source/category counts:
SFR 85, carbon star 27, dark cloud 66, and LOS/diffuse cloud 24. The 2026/current
preview contains SFR 90, carbon star 28, dark cloud 116, and diffuse cloud 24.
The KDE bandwidth follows the legacy convention (`bw = 0.5`). The plot keeps a
source-label override so the 2021 historical view can display `LOS Cloud`,
while modern views use `Diffuse Cloud`.

The accepted production reproduction uses automatic count-label placement near
the KDE curves rather than hard-coded legacy coordinates. This slightly relaxes
pixel-level reproduction but makes the figure robust as source-category
distributions change. Manual `count_label_overrides` remain available for final
manuscript tuning.

For the 2026 manuscript, the accepted production replacement is a horizontal
box/whisker plus jittered exact-value strip plot. It follows the same visual
grammar as the mass-by-wavelength replacement: boxes show the interquartile
range, whiskers span the 10th-90th percentiles, black lines mark medians, and
points show individual detections. The production plot excludes negative-DU
domain artifacts by default because those values come from applying the
neutral/organic DU expression outside its chemically meaningful domain. The
production row order places dark clouds first, followed by carbon stars, SFRs,
and diffuse clouds, so the 2026 plot foregrounds the largest high-DU source
population before the other high-DU source category.

`test_figures_du_by_source_type.py` verifies the 2021 and 2026 source-category
counts, 2021 historical labels, axis labels, KDE construction, production
box/strip construction, and file writing.

## Completed Relative-DU-By-Source-Type Figure Migration

`astromol.figures` now provides `relative_du_by_source_type_data`,
`plot_relative_du_by_source_type`, and `write_relative_du_by_source_type`,
porting the legacy `rel_du_by_source_type` KDE panel through `CensusView`. It
also provides `plot_relative_du_by_source_type_boxplot` and
`write_relative_du_by_source_type_boxplot` as the preferred 2026 production
replacement. The data layer uses secure, non-isotopologue ISM/CSM first
detections by default. Each molecule is credited at most once per generalized
first-detection source category, and fullerenes are excluded by default.

Relative DU is computed as `du / maxdu` using the legacy formula-domain
convention. The 2021 verification data reproduce the legacy source/category
counts: SFR 85, carbon star 27, dark cloud 66, and LOS/diffuse cloud 24. The
2026/current preview contains SFR 90, carbon star 28, dark cloud 116, and
diffuse cloud 24. The legacy KDE panel keeps a source-label override so the
2021 historical view can display `LOS Cloud`, while modern views use
`Diffuse Cloud`.

For the 2026 manuscript, the accepted production replacement is a horizontal
box/whisker plus jittered exact-value strip plot. It uses the same row order
and visual grammar as the absolute DU-by-source box/strip plot: dark clouds,
carbon stars, SFRs, then diffuse clouds. The production view excludes negative
DU-derived domain artifacts by default and pads the left side of the x-axis so
values at zero remain visible instead of sitting on the y-axis spine.

`test_figures_relative_du_by_source_type.py` verifies the 2021 and 2026
source-category counts, value ranges, 2021 historical labels, axis labels, KDE
construction, production box/strip construction, and file writing.

## Completed Mass-By-Source-Type Figure Migration

`astromol.figures` now provides `mass_by_source_type_data`,
`plot_mass_by_source_type`, and `write_mass_by_source_type`, porting the
legacy `mass_by_source_type` KDE figure through `CensusView`. It also provides
`plot_mass_by_source_type_boxplot` and `write_mass_by_source_type_boxplot` as
the preferred 2026 production replacement. The data layer uses secure,
non-isotopologue ISM/CSM first detections by default. Each molecule is credited
at most once per generalized first-detection source category, and fullerenes
are excluded by default so their high masses do not dominate the comparison.

The 2021 verification data reproduce the legacy source/category counts:
SFR 87, carbon star 58, dark cloud 67, and LOS/diffuse cloud 24. The
2026/current preview contains SFR 92, carbon star 71, dark cloud 117, and
diffuse cloud 24. The legacy KDE bandwidth follows the same convention as the
other migrated KDE figures (`bw = 0.5`) and keeps a source-label override so
the historical 2021 view can display `LOS Cloud`, while modern views use
`Diffuse Cloud`.

For the 2026 manuscript, the accepted production replacement is a horizontal
box/whisker plus jittered exact-value strip plot. It uses the same row order
and visual grammar as the DU/source box/strip plots: dark clouds, carbon
stars, SFRs, then diffuse clouds. Boxes show the interquartile range, whiskers
span the 10th-90th percentiles, black lines mark medians, and points show
individual detections.

`test_figures_mass_by_source_type.py` verifies the 2021 and 2026
source-category counts, mass ranges, 2021 historical labels, axis labels, KDE
construction, production box/strip construction, and file writing.

## Completed Wavelength-By-Source-Type Figure Migration

`astromol.figures` now provides `wavelength_by_source_type_data`,
`plot_wavelength_by_source_type`, and `write_wavelength_by_source_type`,
porting the legacy `waves_by_source_type` pie-grid figure through
`CensusView`. It also provides `plot_wavelength_by_source_type_stacked_bar`
and `write_wavelength_by_source_type_stacked_bar` as the preferred 2026
production replacement. The data layer uses secure, non-isotopologue ISM/CSM
first detections by default. Each molecule is credited at most once per
generalized first-detection source category, and every wavelength listed for
that credited first detection is counted. This preserves the legacy convention
that multi-wavelength first detections contribute to multiple wavelength
wedges.

The 2021 verification data reproduce the legacy wedge percentages exactly:
carbon stars are 16.7% cm, 69.7% mm, 1.5% sub-mm, and 12.1% IR; dark clouds
are 67.9% cm and 32.1% mm; LOS/diffuse clouds are 15.4% cm, 7.7% mm, 23.1%
sub-mm, 23.1% IR, and 30.8% UV/Vis; and SFRs are 32.3% cm, 62.5% mm, and
5.2% sub-mm. The 2026/current preview keeps the same pie-grid structure while
using `Diffuse Clouds` for the modern source label.

For the 2026 manuscript, the accepted production replacement is a normalized
horizontal stacked bar chart. It preserves the same wavelength-credit data as
the pie grid but makes the cross-source comparison easier to read, labels the
wavelength-credit denominator for each source category, and omits segment
percentage labels by default because the x-axis already encodes the
proportions.

`test_figures_wavelength_by_source_type.py` verifies the 2021 and 2026
wavelength/source count matrices, 2021 historical labels and percentages,
legend labels, pie-grid construction, production stacked-bar construction, and
file writing.

## Completed Kappa Histogram Figure Migration

`astromol.figures` now provides `kappa_histogram_data`,
`plot_kappa_histogram`, and `write_kappa_histogram`, porting the legacy
`kappas` histogram through `CensusView`. The standard data layer uses secure,
non-isotopologue ISM/CSM molecules. Linear rotors are assigned the physically
appropriate limiting value `kappa = -1` in `Molecule.kappa`; nonlinear
molecules require all three rotational constants. Molecules without a usable
kappa remain excluded from the histogram.

The 2021 verification view contains 221 plotted kappa values out of 240 secure
non-isotopologue ISM/CSM molecules; 118 of the plotted values are linear
rotors and 19 molecules remain without kappa. The 2026/current view contains
305 plotted kappa values out of 325 secure non-isotopologue ISM/CSM molecules;
148 of the plotted values are linear rotors and 20 molecules remain without
kappa.

Category-summary bars, stacked summaries, rug/inset histograms, and KDE
variants were previewed and rejected for the 2026 manuscript because they
either discarded too much of the distribution or smoothed the real `kappa = -1`
point mass into an artificial continuous shoulder. The accepted production
view keeps the histogram representation, uses the standard manuscript aspect
ratio and full boxed axes, and adds a top guide from prolate through
asymmetric to oblate with arrowheads anchored at `kappa = -1` and `kappa = +1`.

`test_figures_kappas.py` verifies the 2021 and 2026 kappa counts, histogram
counts, linear-rotor handling, guide labels, boxed axes, and file writing.
