# AGENTS.md

## What This Project Does

`astromol` is a Python package plus curated scientific database for astronomical
molecule detections. It stores production data in version-controlled JSON and
BibTeX files, loads them into typed models, applies census-aware filtering, and
generates publication products:

- LaTeX scalar fragments and manuscript tables
- Manuscript figures (legacy-repro + modernized production views)
- PowerPoint slides (ISM/CSM and PPD inventories)
- Standard generated-output bundles for GitHub Pages

Primary audience:

- Scientific users who want current figures/slides quickly
- Curators who add/update molecules, detections, sources, telescopes
- Maintainer workflows for reproducible census paper production

---

## Goal Of The Current Refactor (And Why)

The legacy code was one-off and difficult to evolve. This refactor separates:

1. Curated data (`astromol/data/*.json`, `references.bib`)
2. Typed model + validation layer (`models.py`, `database.py`, `validation.py`)
3. Census-view selection logic (`census.py`)
4. Output generation (tables/figures/slides) from shared view logic

Why:

- Make every output reproducible from one authoritative database
- Keep historical census reproduction possible (`2018`, `2021`) while enabling
  dynamic current/`2026` products
- Support regular in-between-census updates without manual paper-era scripts
- Prepare for stable package/API usage and public distribution

---

## Current Architecture (Key Files And Responsibilities)

### Core Data + Model Layer

- `astromol/data/molecules.json`: molecule records (chemistry, properties, tags, manuscript prose, history)
- `astromol/data/detections.json`: detection records (context, source, telescope, status, refs, relationships, history)
- `astromol/data/sources.json`: source objects and classification metadata
- `astromol/data/telescopes.json`: facility objects and metadata
- `astromol/data/references.bib`: Zotero-exported bibliography (not hand-edited)

- `astromol/models.py`: dataclasses + allowed-value constants + derived chemistry properties
- `astromol/database.py`: loads JSON/BibTeX and resolves cross-references with stable IDs
- `astromol/validation.py`: CLI/data integrity checks (`astromol-validate`)

### View/Selection Layer

- `astromol/census.py`:
  - `CensusView.for_census(db, "YYYY")` for frozen census boundaries
  - `CensusView.current(db)` for live/latest boundary
  - Filters for `include_tentative`, `include_disputed`, `include_isotopologues`
  - `FilteredCensusView` for user-defined subset analyses

### Output Generation Layer

- `astromol/latex.py`: scalar fragments + manuscript tables
- `astromol/figures/__init__.py`: public figures API re-export
- `astromol/figures/_core.py`: current figure implementation body (large; still being split)
- `astromol/figures/style.py`: shared style constants (typography, sizing, palette anchors)
- `astromol/slides.py`: slide layout planning, reports, and PowerPoint rendering

### New Canonical Output Inventory Layer

- `astromol/registry.py`:
  - `FIGURE_OUTPUTS`, `TABLE_OUTPUTS`, `SLIDE_OUTPUTS`
  - Stable names, labels, descriptions, and callables for standard outputs
  - `OutputContext` for selected view + baseline view + cached intermediate data

- `astromol/outputs.py`:
  - `astromol-generate-outputs` CLI for standard bundle generation
  - Uses `astromol.registry` as source of truth
  - Writes figures/tables/slides, zip bundle, and public-facing GitHub Pages `index.html`
- `astromol/cli.py`:
  - Public `astromol` CLI for listing and selectively generating standard outputs
  - Delegates complete bundle generation to `astromol.outputs`

### Curation + Docs + Tests

- `scripts/stage_records.py`: staging YAML -> preview/apply workflow with history normalization
- `curation/templates/*.yaml`: full-field curator templates
- `curation/README.md`: staging instructions and quoting rules
- `docs/`: Sphinx + Read the Docs sources
- `docs/notebooks/*.ipynb`: Colab workflows (quickstart, figures, custom views, tables/slides)
- `docs/calculations/*.ipynb`: tracked calculation notebooks used to document
  project-computed scientific values
- `tests/`: regression scripts + focused unit/integration tests
- `.github/workflows/`: CI, package checks, generated-output publishing

---

## Target Architecture (When This Refactor Is Complete)

1. **Stable public API and CLI**
   - Keep module imports stable for users
   - Keep first-pass public CLI aligned with the registry as output products evolve

2. **Output inventory is single-source**
   - `astromol.registry` remains canonical for docs/notebooks/CLI/workflows
   - No duplicated hard-coded output lists

3. **Further figure modularization**
   - Keep `astromol.figures` public surface stable
   - Continue moving families out of `figures/_core.py` into focused modules

4. **Curation schema drift reduced**
   - Templates, staging script, model schema, validation, and docs stay synchronized

5. **Release-grade package polish**
   - PyPI release flow, docs alignment, and branch/default transition cleanup

---

## Patterns And Conventions In This Codebase

### Data/Schema

- Stable ID prefixes are required and must remain deterministic:
  - Molecules: `mol:...`
  - Detections: `det:...`
  - Sources: `src:...`
  - Telescopes: `tel:...`
- Detection statuses are exactly: `secure`, `tentative`, `disputed`
- History event kinds are exactly: `added`, `updated`, `corrected`
- Use `history.introduced` for tracking entry, `history.accepted` for confirmed census membership

### Citation/Reference

- Internal citekeys: `LastName:Year:FirstPage` (suffix only when needed)
- `references.bib` is exported from Zotero and committed as produced

### View/Filtering

- Output logic must consume `CensusView`/`FilteredCensusView`, not ad hoc JSON filters
- Secure/non-isotopologue defaults are conservative; explicit opt-in for tentatives/disputed/isotopologues
- PPD products intentionally include isotopologues by default

### Output Design

- Data builder + render/writer separation is preferred
- Preserve legacy reproducibility modes where intentionally supported
- Modern production style:
  - color-blind-friendlier multi-series palettes
  - shared typography/sizing conventions
  - black axes/full box unless figure-specific reason says otherwise
  - brand anchors: `ASTROMOL_BLUE` default; `MIT_RED` emphasis only
- Generated-output landing page conventions:
  - keep language user-facing rather than implementation-facing
  - surface bundle, key figure, and standard slide decks first
  - internal diagnostics and LaTeX fragments do not need to be public-facing
  - direct image/PDF figure links may open in a new tab/window for convenience

### Testing/Verification

- Keep regression scripts for migration fidelity
- Add focused unit tests for new behavior and refactors
- Validate docs builds (`sphinx -W`) and package builds before major transitions

---

## Decisions Already Made (Do Not Revisit Without Explicit Request)

1. Data is human-curated and human-verified; AI assists code/workflow only.
2. Use MIT license for this project.
3. Base install includes generation deps (matplotlib/scipy/python-pptx), not split into end-user extras.
4. `CensusView` is the central selection mechanism for outputs.
5. Use `history.accepted` to distinguish accepted census membership from discussed-but-unconfirmed records.
6. Keep tentative/disputed detections as first-class detection records.
7. Registry (`astromol.registry`) is canonical for standard output inventories.
8. Keep both historical reproduction capability and modern production outputs.
9. Keep “Diffuse Cloud” terminology (replacing legacy LOS-cloud naming).
10. `references.bib` Zotero export diffs are acceptable and should be committed with relevant work.

---

## Things That Must Not Break

1. `Database()` load and cross-reference resolution across all production JSON and BibTeX data.
2. Census boundary semantics:
   - `for_census("2021")` reproducibility
   - `for_census("2026")` / `current()` behavior
3. Stable detection IDs and reciprocal detection relationship fields.
4. Standard output bundle generation (`astromol-generate-outputs`) and published key assets.
5. PPD slide isotopologue inclusion default.
6. Curation staging workflow (`scripts/stage_records.py`) and template compatibility.
7. Docs + Colab pathways users rely on for figure/slide reproduction.

---

## Where We Are In The Refactor Right Now

Current refactor status:

- Large-scale code changes are **paused for now**:
  - the current priority is scientific database population and cleanup for the
    2026 census rather than further broad architecture work
  - inherited dipole-moment placeholders have been resolved; new curation
    should continue preserving explicit provenance for computed and literature
    values
  - new molecules, detections, source updates, and telescope updates should be
    staged through the existing YAML workflow before production JSON is edited
  - it is acceptable for the maintainer to provide curated data in chat for
    Codex to draft staging YAML, but the maintainer must review and approve the
    staged records before they are applied

- Figure API split is **in progress**:
  - `astromol.figures` moved from single file to package re-export
  - implementation currently in `figures/_core.py`
  - shared style constants moved to `figures/style.py`

- Output registry is **implemented**:
  - `astromol.registry` added
  - `astromol.outputs` now consumes registry specs
  - docs + notebooks now reference registry inventories
  - new tests: `tests/test_registry.py`, `tests/test_outputs.py`

- Generated-output GitHub Pages landing page is **implemented/polished**:
  - `astromol.outputs` now writes a public-facing `index.html` with:
    - a primary-download row for the bundle, cumulative detections figure, and
      the two standard slide decks
    - a multi-column figure-card section below with grouped PNG/PDF links
    - no public-facing LaTeX table links or slide-layout report links
  - PNG and PDF figure links open in a new tab/window
  - page copy was simplified to avoid internal implementation terminology

- First-pass public CLI for standard figures/tables/slides is **implemented**:
  - installed `astromol` command lists registry outputs
  - `astromol figure`, `astromol table`, and `astromol slide` generate selected outputs by stable name
  - `astromol outputs` delegates to the full standard bundle generator
  - committed and pushed as `94cc87e` (`Add public registry-backed CLI`)
  - GitHub Actions for that push completed successfully:
    - `CI`
    - `Package`
    - `Generated Outputs`

- Immediate cleanup items from `CODEBASE_REVIEW.md` are **cleared locally**:
  - source/telescope templates expose `history` and are covered by
    `tests/test_curation_templates.py`
  - `Database()` rejects duplicate source and telescope nicks, covered by
    `tests/test_database.py`
  - `DipoleMoment.total` returns `None` for nonnumeric placeholder components,
    covered by `tests/test_models.py`
  - horizontal boxplots now use Matplotlib `orientation="horizontal"` instead
    of deprecated `vert=False`
  - focused verification passed:
    `python -m pytest tests/test_models.py tests/test_database.py tests/test_curation_templates.py`
  - boxplot regression scripts passed with
    `python -W error::DeprecationWarning`

- Known scientific-data warning backlog:
  - inherited dipole placeholder values (`*`) have been resolved locally;
    validation should report zero warnings unless new curation follow-ups are
    intentionally accepted.

- Working tree was clean after the public CLI work was committed and pushed.
  Current local edits after cleanup are intentionally uncommitted unless the
  user explicitly requests a commit/push.

Refer to [CODEBASE_REVIEW.md](/Users/brett/Dropbox/Programs/census_scripts/astromol/CODEBASE_REVIEW.md) for detailed audit status and verification history.

---

## Rules For Future Work On This Codebase

### Session Start Protocol

At the start of each session, read in this order:

1. `AGENTS.md` (this file)
2. `CODEBASE_REVIEW.md`
3. `README.md` quick status + usage assumptions
4. `SPEC.md` for schema/output semantics
5. `MANUSCRIPT_NOTES_2026.md` for manuscript-impact reminders
6. `git status --short` for in-progress local state

### Implementation Rules

1. Do not invent parallel schema paths. Update model/staging/template/docs together when schema changes.
2. Route new standard outputs through `astromol.registry` and use stable `name` keys.
3. Keep public import paths stable; prefer internal reorganization behind re-exports.
4. Add/adjust tests with every nontrivial refactor, especially for census boundaries and output generation.
5. Keep docs/notebooks in sync when adding/removing standard outputs.

### Response Rules

1. When suggesting a next step for any project task, include the reasoning level
   using exactly: `Recommended Reasoning: X`.

### Curation Rules

1. Stage all new records through YAML templates before applying.
2. Use full templates (all fields visible) for curator ergonomics.
3. Ensure `history` and status semantics are explicit and correct.
4. For molecule promotion changes, prefer one coherent commit per molecule theme.
5. After a staged molecule is promoted and verified, remove obsolete staging YAML unless explicitly preserved for audit.
6. For chat-assisted curation, draft staging YAML only from information the
   maintainer supplies or from explicitly requested source lookups; do not
   invent scientific values, classifications, references, or interpretations.
7. Treat maintainer review as required before running `stage_records.py --apply`
   on newly drafted scientific data.
8. When a database value is computed as part of this work, keep the supporting
   notebook under `docs/calculations/`, cite the software/method references in
   structured `refs` fields, and record the notebook path in a note or history
   summary. Do not track `.ipynb_checkpoints/`.
9. For project-computed dipole moments, add the supporting citekeys to both
   `dipole.refs` and `refs.computation`, record the notebook path in
   `dipole.note`, and use this exact short manuscript sentence in `latex_body`:
   `The dipole moment was calculated as part of this work (\ref{sec:dipole}).`
   Do not repeat the full method/package/basis-set sentence in every molecule;
   that detail belongs in the manuscript dipole section.
10. After any accepted production data change, run
   `python scripts/update_data_baseline.py`, inspect
   `tests/baselines/production_data.json`, and include the intentional baseline
   diff with the curation change.

### Git/Workflow Rules

1. Do not commit/push unless explicitly requested in the current session.
2. Never revert unrelated user changes.
3. Do not use destructive git commands unless explicitly requested.
4. Keep temporary previews/build artifacts out of commits.
5. When producing audit artifacts, document purpose and cleanup path in `SPEC.md`/`CODEBASE_REVIEW.md`.

### Verification Rules

Before declaring substantial work done:

1. Run focused tests for touched area.
2. For production data changes, run
   `python -m pytest tests/test_load.py tests/test_validation.py` after
   refreshing the committed data baseline.
3. Run full `pytest` for broad refactors.
4. Run docs build (`sphinx -W`) when docs/API surface changes.
5. Run clean package build when module/package structure changes.

---

## Useful Commands

```bash
# Core checks
python -m pytest
python -m astromol.validation

# Standard output bundle
astromol-generate-outputs --output-dir build/astromol_outputs --view current --formats png pdf
astromol outputs --output-dir build/astromol_outputs --view current --formats png pdf

# Selective standard outputs
astromol list figures
astromol figure cumulative_detections --view current --output cumulative_detections.pdf
astromol table ism_tables --view 2026 --output-dir build/tables
astromol slide ism_molecule_slide --view current --output-dir build/slides --report

# Docs
python -m sphinx -W -b html docs docs/_build/html

# Curation staging
python scripts/stage_records.py --staging curation/staging/example.yaml
python scripts/stage_records.py --staging curation/staging/example.yaml --apply
python scripts/update_data_baseline.py
```
