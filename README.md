# astromol

`astromol` is a Python package and curated database of astronomical molecule
detections. You can use it to query the inventory and generate census figures,
tables, and PowerPoint slides from the same records. Molecules, detections,
sources, telescopes, and references are stored in version-controlled files;
the Python API loads them and checks the links between records.

## Latest Figures And Slides

The latest standard figures and PowerPoint slides are generated automatically
from the current database by GitHub Actions and published to GitHub Pages:

<https://bmcguir2.github.io/astromol/>

If you need the latest figure or slide, start there. Direct downloads are also
linked below:

- [latest ISM/CSM cumulative detections figure](https://bmcguir2.github.io/astromol/figures/png/cumulative_detections.png)
- [latest ISM/CSM detections slide](https://bmcguir2.github.io/astromol/slides/astro_molecules_current.pptx)
- [latest protoplanetary disk detections slide](https://bmcguir2.github.io/astromol/slides/ppd_molecules_current.pptx)
- [complete latest output bundle](https://bmcguir2.github.io/astromol/astromol_latest_outputs.zip)

To choose a different census view, filter the database, or regenerate individual
outputs, use the Google Colab notebooks linked in the documentation. Outputs
from pull requests and branch builds are also available from the
[Generated Outputs workflow](https://github.com/bmcguir2/astromol/actions/workflows/generated-outputs.yml).

This README describes the active `refactor` branch. The database, curation
workflow, and output tools are usable. The package and public API are still
pre-release, and this branch has not yet replaced `main`.

## Current Status

- Production data lives in `astromol/data/*.json`.
- Bibliographic metadata lives in `astromol/data/references.bib`.
- `references.bib` is exported from Zotero and should not be hand-edited.
- `SPEC.md` is the source of truth for the current schema and architecture.
- `MANUSCRIPT_NOTES_2026.md` collects writing-time reminders for the next
  census paper (originally planned for 2026).
- New records should be staged through YAML files in `curation/staging/` before
  being applied to production JSON.
- Work is now focused on preparing the database for the next census. Stage new
  molecule and detection updates for maintainer review before applying them.
  Broad refactor work is paused.
- Calculation notebooks that support project-computed database values are
  tracked under `docs/calculations/`; local `.ipynb_checkpoints/` directories
  are ignored.

## Quick Load Check

```bash
python -m pytest tests/test_load.py
```

Or from Python:

```python
from astromol.database import Database

db = Database()
print(len(db.molecules), len(db.detections))
```

`Database()` loads the five data files, links detections to molecule, source,
telescope, and reference records, and checks detection IDs.

## Census Views And Output Generation

Use a `CensusView` to choose which records belong in an output. The same
membership and filtering rules then apply to figures, tables, slides, and
manuscript counts:

```python
from pathlib import Path

from astromol.census import CensusView
from astromol.database import Database
from astromol.figures import (
    cumulative_detection_data,
    write_cumulative_detections_plot,
)
from astromol.slides import write_molecule_slide, write_ppd_detection_slide

db = Database()
view = CensusView.current(db)
data = cumulative_detection_data(view)
write_cumulative_detections_plot(data, Path("cumulative_detections.pdf"))
write_molecule_slide(view, Path("astro_molecules.pptx"), profile="balanced")
write_ppd_detection_slide(view, Path("ppd_molecules.pptx"))
```

Use `CensusView.for_census(db, "2018")` or `"2021"` for membership in a
published census, and `CensusView.current(db)` for the working inventory.
The next census is still in preparation; its eventual publication year does
not define which records belong in it. `for_census(db, "2026")` remains an
alias for `current` for compatibility. Historical selection restores membership;
exact reproduction of an old output requires its original revision and
dependency environment.

Figure helpers live in `astromol.figures`. They separate data selection from
plotting so the selected records and calculated values can be checked before
rendering. PowerPoint helpers live in `astromol.slides`; the ISM/CSM slide has
a legacy 2021 layout and a balanced layout that adjusts to the inventory and
reports layout problems. See `SPEC.md` for the output API and
`GENERATION_MIGRATION.md` for the checks against the 2021 census.
The PPD slide helper includes detected isotopologues by default.

Generated slides distinguish the software version from the database freshness
date. The top-right credit line reports the installed package version when
available and includes a traceable git hash when run from a checkout. The lower
count/date block reports the latest modification date across the selected
molecules, detections, sources, and telescopes.

The standard latest output bundle can be regenerated locally with:

```bash
astromol-generate-outputs --output-dir build/astromol_outputs --view current --formats png pdf
```

To list available outputs or generate one figure, table group, or slide deck,
use the `astromol` command:

```bash
astromol list figures
astromol figure cumulative_detections --view current --output cumulative_detections.pdf
astromol table ism_tables --view current --output-dir build/tables
astromol slide ism_molecule_slide --view current --output-dir build/slides --report
astromol outputs --output-dir build/astromol_outputs --view current --formats png pdf
```

## Curation Workflow

Copy a template from `curation/templates/` into `curation/staging/`, fill in the
fields you know, then generate a preview:

```bash
python scripts/stage_records.py --staging curation/staging/example.yaml
```

The report's `Generated Count Updates` section shows how the batch would change
the inventory and regression baseline. Review the report, full preview JSON,
and proposed baseline under `astromol/data/`. Once the maintainer has approved
them, apply the records:

```bash
python scripts/stage_records.py --staging curation/staging/example.yaml --apply
```

Applying staged records refreshes
`tests/baselines/production_data.json` automatically. Do not manually edit
generated count expectations in regression scripts; update the curated data and
let the staging workflow regenerate the baseline.

To update an existing record, first generate its full-field, digest-locked
template:

```bash
python scripts/stage_records.py \
  --prepare-update molecule mol:EXAMPLE \
  --output curation/staging/example_update.yaml
```

Additions and updates can share a batch. Declare a detection relationship on
the scientifically meaningful side; the script adds the reciprocal link and
shows it in the preview. Changes to scientific status or census acceptance
must be explicit in the staged records and reviewed by the curator. The
complete preview must pass validation before apply.

Successful staging runs also write
`astromol/data/<name>_stage_manifest.json`. After the curation batch is
finished, use the cleanup helper to remove the staging YAML, preview JSON,
stage report, and manifest:

```bash
python scripts/cleanup_stage.py --name example
```

If you want the cleanup helper to also create a git commit, it automatically
stages the production JSON files touched by the manifest plus modified
`astromol/data/references.bib` and
`tests/baselines/production_data.json` by default:

```bash
python scripts/cleanup_stage.py --name example --commit-message "Add Example curation batch"
```

Commit-mode cleanup verifies the applied manifest hashes and runs
`python scripts/check_curation.py` before deleting review artifacts or
committing.

If the curation batch resolves one or more tracked GitHub issues, pass
`--close-issue` with `--push`. The helper closes each issue after the push
succeeds and comments with a link to the new commit:

```bash
python scripts/cleanup_stage.py \
  --name example \
  --commit-message "Add Example curation batch" \
  --push \
  --close-issue 123
```

After an approved production data change, inspect the generated baseline diff
before committing, then run the local curation verification check:

```bash
python scripts/check_curation.py
```

This runs production-data validation plus the curation-sensitive load,
validation-baseline, and output regression tests. The output regression step is
needed because new secure molecules/detections can legitimately change current
scientific counts and the space required by tables, figures, and slides.

See `curation/README.md` for template notes, YAML quoting rules, and schema
maintenance expectations.

For assisted curation, draft YAML from curator-provided data or explicitly
requested source lookups. The maintainer must check the scientific fields,
history, and references before apply. If any reviewed input changes, generate
and review a fresh preview.

When a database value is calculated as part of this project, keep the
calculation notebook under `docs/calculations/`, cite the relevant software,
method, and basis-set references in the molecule metadata, and record the
notebook path in a note or history entry for provenance. For project-computed
dipole moments, molecule prose should use the short manuscript note
`The dipole moment was calculated as part of this work (\ref{sec:dipole}).`
The detailed method, package, and basis-set description belongs in the
manuscript dipole-methods section.

## AI Assistance and Data Curation

This project uses OpenAI Codex as a coding assistant for software
implementation, refactoring, validation scripts, documentation drafts, and
workflow support.

Scientific database content is not AI-generated. Records and curation decisions
come from the literature, Zotero-managed bibliography exports, legacy census
materials, and human review.
AI assistance can help draft staging files or transform supplied records, but
the maintainer must verify all scientific data before they enter production.

The project maintainer is responsible for all committed code, data, and
documentation.

## Reference Workflow

References are maintained in Zotero and exported to
`astromol/data/references.bib`. When new data records use new citekeys, export
the updated Zotero collection and commit the resulting BibTeX file.

Internal citekeys use:

```text
LastName:Year:FirstPage
```

Add lowercase suffixes only when needed to disambiguate duplicate keys.

## Important Files

- `SPEC.md`: architecture and schema specification
- `MANUSCRIPT_NOTES_2026.md`: 2026 manuscript reminders and explanatory notes
- `GENERATION_MIGRATION.md`: output-generation migration and verification log
- `astromol/cli.py`: public command-line interface for standard outputs
- `astromol/models.py`: dataclasses and allowed values
- `astromol/database.py`: JSON/BibTeX loader and cross-reference resolver
- `astromol/figures/`: figure data builders and plotting helpers
- `astromol/registry.py`: canonical registry of standard figures, tables, and slides
- `astromol/outputs.py`: standard generated-output bundle writer
- `astromol/slides.py`: PowerPoint slide builders and layout diagnostics
- `astromol/data/`: production data files
- `curation/templates/`: curator-facing YAML templates
- `docs/`: Sphinx/Read the Docs documentation source
- `docs/calculations/`: tracked notebooks for project-computed scientific values
- `scripts/stage_records.py`: staging, validation, preview, and apply workflow
- `scripts/update_data_baseline.py`: manually refreshes the generated
  production-data, regression-count, and validation-warning baseline
- `tests/baselines/production_data.json`: expected production inventory and
  curation-sensitive regression-count snapshot used by tests

## Packaging

This branch now includes initial Python packaging metadata in `pyproject.toml`.
For local development, install the package in editable mode with the optional
dependencies you need:

```bash
python -m pip install -e ".[dev]"
```

The base install includes the database loader, census/table helpers, figure
generation, PowerPoint slide generation, and bundled production data. The YAML
curation templates and staging script are contributor-facing source-checkout
tools. Optional dependency groups are limited to contributor-oriented tooling:

- `docs`: Sphinx/Read the Docs documentation build dependencies
- `dev`: test, build, and package-check dependencies

The pre-release package version is `2026.0.0.dev0`, which sorts after the legacy
PyPI release (`2021.7`). It identifies the software and API. Census membership
and database update dates are tracked separately in record history.

## Documentation

Documentation source files live under `docs/` and are configured for Sphinx and
Read the Docs. Install the documentation dependencies with:

```bash
python -m pip install -e ".[docs]"
```

Build the HTML docs locally with:

```bash
python -m sphinx -b html docs docs/_build/html
```

The generated local entry point is `docs/_build/html/index.html`.

Example notebooks live under `docs/notebooks/` and are linked from the docs.
They include Google Colab badges so users can run quickstart, figure, custom
view, table, and slide examples directly from GitHub. During the refactor they
install from the `refactor` branch; after PyPI release, their setup cells should
be updated to install `astromol` from PyPI.

## FAQ

These links open the draft refactor-branch notebooks directly in Google Colab.
Set `VIEW_CHOICE = "current"` in the notebook configuration cells when you want
the live database rather than a historical census view.

### How Do I Get The Latest ISM/CSM Cumulative Detections Figure?

Open the
[Reproduce figures notebook](https://colab.research.google.com/github/bmcguir2/astromol/blob/refactor/docs/notebooks/02_reproduce_figures.ipynb),
run the setup and configuration cells, set `VIEW_CHOICE = "current"`, then run
the **Cumulative ISM/CSM detections** cell. The final notebook cell can download
the generated PNG/PDF outputs.

### How Do I Get The Latest ISM/CSM Detections Slide?

Open the
[Tables and slides notebook](https://colab.research.google.com/github/bmcguir2/astromol/blob/refactor/docs/notebooks/04_tables_and_slides.ipynb),
run the setup and configuration cells, set `VIEW_CHOICE = "current"`, then run
the PowerPoint slide cell. It writes `astro_molecules_current.pptx` using the
production balanced layout.

### How Do I Get The Latest Protoplanetary Disk Detections Slide?

Use the same
[Tables and slides notebook](https://colab.research.google.com/github/bmcguir2/astromol/blob/refactor/docs/notebooks/04_tables_and_slides.ipynb)
with `VIEW_CHOICE = "current"`, then run the PowerPoint slide cell. It writes
`ppd_molecules_current.pptx`; detected isotopologues are included by default.

### How Do I Get Any Of The Other Figures?

Open the
[Reproduce figures notebook](https://colab.research.google.com/github/bmcguir2/astromol/blob/refactor/docs/notebooks/02_reproduce_figures.ipynb),
run the setup/configuration cells, then run the specific figure cells you want.
For custom subsets, start from the
[Custom filtered views notebook](https://colab.research.google.com/github/bmcguir2/astromol/blob/refactor/docs/notebooks/03_custom_views.ipynb).

## Verification Scripts

Run the production-data validator before committing data changes:

```bash
astromol-validate
```

or, from a source checkout without installing entry points:

```bash
python -m astromol.validation
```

For production curation changes, run the local curation verification check:

```bash
python scripts/check_curation.py
```

Run the full regression suite with:

```bash
python -m pytest
```

The tests live under `tests/` and verify database loading, census-view
membership, figure/table data builders, slide layout/rendering properties, and
selected 2021/current output regressions. The current validator reports no
accepted production-data warnings.

Pytest also runs the regression scripts used to check the migration. Keeping
those checks preserves the audited comparisons with older outputs. As the API,
examples, and CI workflow stabilize, they can be moved into focused unit and
integration tests, with slow figure and slide checks marked separately.

## Continuous Integration

GitHub Actions is configured in `.github/workflows/ci.yml`. On pushes, pull
requests, and manual workflow dispatches, CI installs the package with
development and documentation dependencies, runs `astromol-validate`, runs the
pytest suite, and builds the Sphinx documentation with warnings treated as
errors.
