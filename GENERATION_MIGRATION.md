# Census Output Generation Migration

This log tracks modernization of the legacy census-output functions from
`astromol/functions_legacy.py` into the refactored database-backed package API.
It is a working audit artifact for the migration; production scientific data
remain the JSON files under `astromol/data/`.

## Reference Inputs

- `astromol/functions_legacy.py`: legacy 2021 output-generation code.
- `astromol/2021_census_arxiv.tex`: 2021 manuscript source used as a reference
  for generated inputs and manuscript integration.

The legacy module is treated as a behavioral specification, not as runnable
production code. It imports removed legacy modules (`astromol.molecules`,
`astromol.sources`, `astromol.telescopes`) and contains many embedded LaTeX
strings that emit Python syntax warnings under modern parsing.

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
- The PowerPoint molecule slide is a first-class product and will be ported.
- A PowerPoint slide summarizing PPD detections should also be added as a
  first-class product.
- The migration will proceed one function group at a time, with parity checks
  against the 2021 behavior where possible.

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
| 4 | Port figure data builders and plots. | Pending |
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
| `cumu_det_plot` | 276-478 | Figure | `astromol.figures` | Pending | Visual comparison to 2021 figure |
| `cumu_det_natoms_plot` | 481-747 | Figure | `astromol.figures` | Pending | Visual comparison to 2021 figure |
| `det_per_year_per_atom` | 750-867 | Figure | `astromol.figures` | Pending | Visual comparison to 2021 figure |
| `facility_shares` | 870-1044 | Figure | `astromol.figures` | Pending | Visual comparison to 2021 figure |
| `scopes_by_year` | 1047-1202 | Figure | `astromol.figures` | Pending | Visual comparison to 2021 figure |
| `periodic_heatmap` | 1205-1436 | Figure | `astromol.figures` | Pending | Visual comparison to 2021 figure |
| `mass_by_wavelengths` | 1439-1673 | Figure | `astromol.figures` | Pending | Visual comparison to 2021 figure |
| `mols_waves_by_atoms` | 1676-1818 | Figure | `astromol.figures` | Pending | Visual comparison to 2021 figure |
| `du_histogram` | 1821-1894 | Figure | `astromol.figures` | Pending | Visual comparison to 2021 figure |
| `type_pie_chart` | 1897-2114 | Figure | `astromol.figures` | Pending | Visual comparison to 2021 figure |
| `source_pie_chart` | 2117-2311 | Figure | `astromol.figures` | Pending | Visual comparison to 2021 figure |
| `indiv_source_pie_chart` | 2314-2500 | Figure | `astromol.figures` | Pending | Visual comparison to 2021 figure |
| `mol_type_by_source_type` | 2503-2697 | Figure | `astromol.figures` | Pending | Visual comparison to 2021 figure |
| `du_by_source_type` | 2700-2850 | Figure | `astromol.figures` | Pending | Visual comparison to 2021 figure |
| `rel_du_by_source_type` | 2853-2992 | Figure | `astromol.figures` | Pending | Visual comparison to 2021 figure |
| `mass_by_source_type` | 2995-3152 | Figure | `astromol.figures` | Pending | Visual comparison to 2021 figure |
| `waves_by_source_type` | 3155-3356 | Figure | `astromol.figures` | Pending | Visual comparison to 2021 figure |
| `kappas` | 3359-3417 | Figure | `astromol.figures` | Pending | Visual comparison if retained |
| `waves_pie_chart` | 3420-3556 | Figure | `astromol.figures` | Pending | Visual comparison if retained |
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
- Re-evaluate `periodictable`, `colour`, `seaborn`, and `scipy` as each plot is
  ported; avoid carrying dependencies that only supported legacy styling.

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

## Immediate Next Step

Begin porting figure generators. The first candidate is `cumu_det_plot`, which
should consume `CensusView`-derived first-detection years and produce both a
2021 verification preview and a current 2026 manuscript preview.
