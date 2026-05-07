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

- Output-generation code should live inside the installable `astromol` API.
- Reusable census-selection logic will be centralized before porting output
  functions. Tables, figures, and slides must consume the same view/filtering
  layer rather than independently filtering raw records.
- Users should eventually be able to generate an up-to-date complete compiled
  PDF, figures, tables, and the summary slide from the current database between
  formal publication releases.
- The PowerPoint molecule slide is a first-class product and will be ported.
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
| 2 | Port scalar LaTeX generators. | Pending |
| 3 | Port LaTeX table generators. | Pending |
| 4 | Port figure data builders and plots. | Pending |
| 5 | Port PowerPoint molecule slide generation. | Pending |
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
| `make_ism_tables` | 3563-3720 | LaTeX table | `astromol.latex` | Pending | Reproduce 2021 ISM tables |
| `make_exgal_count` | 3722-3739 | LaTeX scalar | `astromol.latex` | Pending | Reproduce `nexgal.tex` value |
| `make_exgal_percent` | 3741-3758 | LaTeX scalar | `astromol.latex` | Pending | Reproduce `nexgalpercent.tex` value |
| `_make_exgal_sentence` | 3760-3783 | LaTeX prose fragment | `astromol.latex` | Pending | Reproduce exgal source sentence |
| `make_det_count` | 3785-3802 | LaTeX scalar | `astromol.latex` | Pending | Reproduce `ndetects.tex` value |
| `make_elem_count` | 3805-3832 | LaTeX scalar | `astromol.latex` | Pending | Reproduce `nelems.tex` value |
| `make_ppd_count` | 3834-3859 | LaTeX scalar | `astromol.latex` | Pending | Reproduce `nppds.tex` value |
| `make_exo_count` | 3861-3885 | LaTeX scalar | `astromol.latex` | Pending | Reproduce `nexos.tex` value |
| `make_ices_count` | 3887-3911 | LaTeX scalar | `astromol.latex` | Pending | Reproduce `nices.tex` value |
| `make_ppd_isos_count` | 3913-3938 | LaTeX scalar | `astromol.latex` | Pending | Reproduce `nppdisos.tex` value |
| `make_exgal_table` | 3941-4140 | LaTeX table | `astromol.latex` | Pending | Reproduce 2021 exgal table |
| `make_ppd_table` | 4143-4301 | LaTeX table | `astromol.latex` | Pending | Reproduce 2021 PPD table |
| `make_exo_table` | 4304-4406 | LaTeX table | `astromol.latex` | Pending | Reproduce 2021 exoplanet table |
| `make_ice_table` | 4408-4523 | LaTeX table | `astromol.latex` | Pending | Reproduce 2021 ice table |
| `make_det_per_year_by_atoms_table` | 4525-4655 | LaTeX table | `astromol.latex` | Pending | Reproduce rate table |
| `make_facility_table` | 4657-4753 | LaTeX table | `astromol.latex` | Pending | Reproduce 2021 facility table |
| `make_source_table` | 4756-4860 | LaTeX table | `astromol.latex` | Pending | Reproduce 2021 source table |
| `make_rate_counts` | 4862-4913 | LaTeX scalar group | `astromol.latex` | Pending | Reproduce rate scalar files |
| `make_percent_radio` | 4915-4939 | LaTeX scalar | `astromol.latex` | Pending | Reproduce `radiopercent.tex` value |
| `make_scopes_count` | 4941-4964 | LaTeX scalar | `astromol.latex` | Pending | Reproduce `nscopes.tex` value |
| `make_percent_unsat` | 4966-4995 | LaTeX scalar | `astromol.latex` | Pending | Reproduce `unsatpercent.tex` value |
| `make_sat_list` | 4997-5026 | LaTeX scalar/list | `astromol.latex` | Pending | Reproduce `satlist.tex` content |
| `make_sat_count` | 5028-5054 | LaTeX scalar | `astromol.latex` | Pending | Reproduce `nsats.tex` value |
| `make_sat_percent` | 5056-5085 | LaTeX scalar | `astromol.latex` | Pending | Reproduce `satpercent.tex` value |
| `make_sfr_rad_percent` | 5087-5113 | LaTeX scalar | `astromol.latex` | Pending | Reproduce `sfr_rad_percent.tex` value |
| `make_dark_rad_percent` | 5115-5141 | LaTeX scalar | `astromol.latex` | Pending | Reproduce `dark_rad_percent.tex` value |
| `make_mols_slide` | 5147-5561 | PowerPoint slide | `astromol.slides` | Pending | Visual comparison to legacy slide |

## Dependency Notes

Legacy imports include `numpy`, `matplotlib`, `periodictable`, `colour`,
`seaborn`, `scipy`, and `python-pptx`. The modern implementation should add
project dependency metadata before these become required package features.
Each dependency should be justified by the migrated code:

- Keep `numpy` and `matplotlib` for plotting unless a specific replacement is
  chosen.
- Keep `python-pptx` for the molecule slide product.
- Re-evaluate `periodictable`, `colour`, `seaborn`, and `scipy` as each plot is
  ported; avoid carrying dependencies that only supported legacy styling.

## Completed View Verification

`test_census_view.py` verifies the shared view layer against the audited
historical table memberships:

- 2018 ISM/CSM: 204 secure molecules
- 2018 extragalactic: 63 secure, 65 with tentative rows
- 2018 exoplanet: 5 secure molecules
- 2018 PPD: 35 secure molecules
- 2021 ISM/CSM: 240 secure molecules
- 2021 extragalactic: 73 secure, 75 with tentative rows
- 2021 exoplanet: 9 secure molecules
- 2021 PPD: 40 secure molecules
- 2021 source/facility contribution counts
- 2026 census view and current view currently return identical secure
  accepted detection sets by context

## Immediate Next Step

Begin porting scalar LaTeX generators, starting with `make_det_count`,
`make_elem_count`, and the context count functions. These should consume
`CensusView` rather than filtering raw database records.
