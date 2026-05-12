# astromol

`astromol` is a Python library and curated data set for cataloging astronomical
molecule detections. It stores molecule, detection, source, telescope, and
reference metadata in version-controlled data files, then loads and resolves
those records through a small Python API.

This README describes the active refactor branch. The data model and curation
workflow are usable, but the public package/API documentation is still being
developed before this branch replaces `main`.

## Current Status

- Production data lives in `astromol/data/*.json`.
- Bibliographic metadata lives in `astromol/data/references.bib`.
- `references.bib` is exported from Zotero and should not be hand-edited.
- `SPEC.md` is the source of truth for the current schema and architecture.
- `MANUSCRIPT_NOTES_2026.md` collects writing-time reminders for the 2026
  census paper.
- New records should be staged through YAML files in `curation/staging/` before
  being applied to production JSON.

## Quick Load Check

```bash
python test_load.py
```

Or from Python:

```python
from astromol.database import Database

db = Database()
print(len(db.molecules), len(db.detections))
```

`Database()` loads references, telescopes, sources, molecules, and detections;
it also resolves cross-references and validates stable detection IDs.

## Census Views And Output Generation

Tables, scalar manuscript fragments, and figures should be generated through a
`CensusView`, which applies the accepted census boundary and the standard
isotopologue/tentative/disputed filtering rules consistently:

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
view = CensusView.for_census(db, "2026")
data = cumulative_detection_data(view)
write_cumulative_detections_plot(data, Path("cumulative_detections.pdf"))
write_molecule_slide(view, Path("astro_molecules.pptx"), profile="balanced")
write_ppd_detection_slide(view, Path("ppd_molecules.pptx"))
```

Use `CensusView.for_census(db, "2021")` for historical reproduction and
`CensusView.current(db)` for the live database. Until the 2026 census cutoff is
frozen, the 2026 and current views are expected to match.

Figure helpers live in `astromol.figures`. Each migrated figure has a data
builder, plotting function, and writer function so the scientific selection can
be tested separately from the visual rendering. PowerPoint helpers live in
`astromol.slides`; the current molecule-slide implementation includes both the
legacy 2021 layout profile and a balanced dynamic profile with layout
diagnostics. See `SPEC.md` for the complete list of table, figure, and slide
helpers, and `GENERATION_MIGRATION.md` for the verification/audit trail against
the 2021 census.
The PPD slide helper includes detected isotopologues by default.

Generated slides distinguish the software version from the database freshness
date. The top-right credit line reports the installed package version when
available, or a traceable development git hash when run from the refactor
checkout. The lower count/date block reports the latest modification date from
the selected curated records.

## Curation Workflow

Copy a template from `curation/templates/` into `curation/staging/`, fill in the
fields you know, then generate a preview:

```bash
python scripts/stage_records.py --staging curation/staging/example.yaml
```

After reviewing the generated report and preview JSON under `astromol/data/`,
apply the staged records:

```bash
python scripts/stage_records.py --staging curation/staging/example.yaml --apply
```

See `curation/README.md` for template notes, YAML quoting rules, and schema
maintenance expectations.

## AI Assistance and Data Curation

This project uses OpenAI Codex as a coding assistant for software
implementation, refactoring, validation scripts, documentation drafts, and
workflow support.

Scientific database content is not AI-generated. Molecule records, detection
records, source/telescope metadata, reference mappings, classifications, and
curation decisions are derived from the scientific literature, Zotero-managed
bibliography exports, legacy census materials, and human curator review. AI
assistance may be used to stage or transform records for review, but no staged
scientific data are accepted into production without human verification.

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
- `astromol/models.py`: dataclasses and allowed values
- `astromol/database.py`: JSON/BibTeX loader and cross-reference resolver
- `astromol/figures.py`: figure data builders and plotting helpers
- `astromol/slides.py`: PowerPoint slide builders and layout diagnostics
- `astromol/data/`: production data files
- `curation/templates/`: curator-facing YAML templates
- `docs/`: Sphinx/Read the Docs documentation source
- `scripts/stage_records.py`: staging, validation, preview, and apply workflow

## Packaging

This branch now includes initial Python packaging metadata in `pyproject.toml`.
For local development, install the package in editable mode with the optional
dependencies you need:

```bash
python -m pip install -e ".[dev]"
```

The base install includes the database loader, census/table helpers, figure
generation, PowerPoint slide generation, and YAML curation workflow. Optional
dependency groups are limited to contributor-oriented tooling:

- `docs`: Sphinx/Read the Docs documentation build dependencies
- `dev`: test-runner dependencies

The package version is currently `0.1.0.dev0` while the refactor branch is still
pre-release. Census boundaries and database freshness are tracked separately in
the record history metadata; the Python package version is the software/API
version.

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

Run the regression suite with:

```bash
python -m pytest
```

The tests live under `tests/` and verify database loading, census-view
membership, figure/table data builders, slide layout/rendering properties, and
selected 2021/current output regressions. The current validator reports known
legacy `*` dipole placeholders as warnings until those values are reconciled.

The current pytest suite intentionally includes a wrapped regression-script
harness that preserves the migration checks used during the refactor. Once the
public API, documentation examples, and CI workflow stabilize, those checks
should be incrementally converted into a more conventional unit/integration
test suite with fixtures, smaller focused assertions, and explicit slow-test
markers for figure and slide generation.

## Continuous Integration

GitHub Actions is configured in `.github/workflows/ci.yml`. On pushes, pull
requests, and manual workflow dispatches, CI installs the package with
development and documentation dependencies, runs `astromol-validate`, runs the
pytest suite, and builds the Sphinx documentation with warnings treated as
errors.
