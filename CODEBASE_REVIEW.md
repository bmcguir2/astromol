# astromol Codebase Review

Date: 2026-05-13
Baseline commit: `6783be5` (`Polish package release metadata`)

Update: large-scale refactor and release-polish work is paused for now. The
active project phase is 2026 scientific database curation: first resolving the
remaining dipole-moment placeholders, then staging new molecules and/or
detections through the existing YAML workflow before applying them to
production JSON.

This review is a fresh pass over the refactor branch after the data model,
curation workflow, table/figure/slide generation, packaging metadata,
documentation, Read the Docs, Colab notebooks, and CI were brought online.

## Executive Summary

The project is in much better shape than the original legacy workflow. The
core split is right: curated JSON data are loaded into typed records,
`CensusView` owns historical/current selection rules, and tables, figures, and
slides consume those views rather than each reimplementing census logic. The
database is loadable from an installed wheel, semantic validation passes, and
the output-generation regression suite is broad enough to catch many accidental
behavior changes.

The main remaining risk for the next phase is curation drift: new database
updates must keep model schema, staging templates, validation, references, and
history semantics aligned. The current codebase is good enough for the next
round of data population, so broad architecture work should wait unless a
curation task exposes a concrete blocker.

## What Works Well

- **Data/view separation is sound.** `Database` resolves records and
  references; `CensusView` handles accepted/current/census logic,
  tentative/disputed inclusion, and isotopologue filtering.
- **Stable IDs and history metadata are now doing useful work.** Detection IDs,
  reciprocal detection relationships, and `history.accepted`/`introduced`
  enable distinctions that the 2021 code could not express cleanly.
- **Output generators are scientifically testable.** The data-builder /
  renderer / writer pattern in `astromol.figures`, `astromol.latex`, and
  `astromol.slides` is the correct direction.
- **Slide generation is more auditable than the legacy code.** Layout planning,
  reports, warnings, and rendering are separated well enough to debug crowded
  layouts before writing binary PowerPoint files.
- **Documentation and examples are no longer afterthoughts.** Read the Docs,
  Colab notebooks, FAQ, and local recipes are in place.
- **The current validation/test baseline is meaningful.** `astromol-validate`
  reports 0 errors and 3 known dipole-placeholder warnings; the pytest
  regression harness exercises table, figure, and slide outputs.

## Priority Findings

### 1. Curation Schema Drift Is The Biggest Practical Risk

The curation schema is currently repeated in several places:

- dataclasses and allowed values in `astromol/models.py`
- staging field order/defaults/required fields in `scripts/stage_records.py`
- templates in `curation/templates/*.yaml`
- prose schema in `SPEC.md` and docs
- validation in `astromol/validation.py`

Previously observed source/telescope template drift has been corrected:

- `curation/templates/source.yaml` and `curation/templates/telescope.yaml`
  expose `history` blocks that match `scripts/stage_records.py`.
- `curation/templates/telescope.yaml` no longer refers to a nonexistent
  `TELESCOPE_TYPES` constant.
- `tests/test_curation_templates.py` checks template keys against
  `stage_records.py` `FIELDS`, `DEFAULTS`, and `REQUIRED`.

Recommended direction for the current curation phase: keep the hand-written
templates and existing consistency tests. Do not introduce generated templates
or a new schema-abstraction layer unless repeated curation work exposes a
specific blocker. Avoid schema drift by updating templates, staging defaults,
validation, specs, and focused template tests in the same change whenever
curator-facing fields change.

### 2. Source And Telescope Duplicate Keys Are Checked

`Database._load_molecules`, `_load_detections`, `_load_sources`, and
`_load_telescopes` now reject duplicate primary keys before assigning records
into lookup dictionaries. Duplicate `source.nick` and `telescope.nick` values
raise `ValueError` during `Database()` load.

Focused coverage: `tests/test_database.py` exercises duplicate source and
telescope nick rejection against a minimal temporary data directory.

### 3. Known Dipole Placeholders Are Warnings And API-Safe

The validator correctly reports the remaining inherited `*` dipole placeholders as
warnings. `DipoleMoment.total` now returns `None` when any populated component
is nonnumeric, so these warning records no longer raise a raw `TypeError` when
users access the total dipole moment.

Current warning records:

- `mol:SO+`
- `mol:MgCN`
- `mol:HNCS`

Recommended action: keep the manuscript-note reminder to resolve the values
before final dipole-based analysis or manuscript claims.

For values computed as part of this work, keep the supporting notebooks under
`docs/calculations/`, cite method/software references in structured molecule
reference fields, and record the local notebook path in the dipole note or
history summary. `docs/calculations/.ipynb_checkpoints/` should not be tracked.

### 4. `astromol.figures` Is Working But Too Large To Scale Comfortably

`astromol/figures.py` is about 7,500 lines, with data classes, style constants,
aggregation logic, plot rendering, and file writers all in one module. This is
manageable right now because the migration is fresh, but it will become a pain
point as new figures and 2026 manuscript variants accumulate.

Recommended direction:

- keep the public import surface stable by re-exporting from `astromol.figures`;
- split implementation into focused modules such as:
  - `figures/style.py`
  - `figures/cumulative.py`
  - `figures/chemistry.py`
  - `figures/source_type.py`
  - `figures/wavelength.py`
  - `figures/facility.py`
  - `figures/registry.py`
- add a central figure registry that notebooks/docs can query instead of
  manually maintaining lists of available figures.

This is not an emergency, but it should happen before many new figures are
added.

### 5. Filtered Views Duplicate Too Much `CensusView` Logic

`FilteredCensusView` wraps `CensusView`, which is the right idea, but it repeats
many convenience methods and count methods. That duplication is currently
small, but every new view method will need to be added twice.

Recommended action: factor shared count/context helpers so both view classes
reuse the same implementation, or make `FilteredCensusView` delegate unknown
attributes to `base` and override only filtering-specific methods.

### 6. Regression Coverage Is Broad But Not Yet A Polished Test Suite

The current tests preserve migration checks, which was the right immediate
choice. The next testing step should convert high-value checks into focused
pytest unit/integration tests:

- formula parsing and isotope mass fallback;
- charge/radical/DU/kappa computations;
- `CensusView` acceptance boundaries, tentative/disputed handling, and
  isotopologue filtering;
- duplicate source/telescope/molecule/detection keys;
- staging defaults/history normalization/reference checking;
- slide layout warnings for crowded inventories.

Also add explicit slow-test markers for full figure/slide generation so local
development can run a fast subset when needed.

### 7. CI Should Add Packaging Verification Before Release

Current GitHub Actions validates data, runs pytest, and builds docs. Before a
public package release, CI should also:

- build sdist and wheel;
- run `twine check`;
- install the wheel into a clean environment;
- load `Database()`;
- run `astromol-validate` from the installed wheel.

That mirrors the manual package audit already performed and protects against
missing package data or import-only-from-source mistakes.

### 8. Public CLI Has A First-Pass Output Surface

The package now exposes both the existing focused commands and a public
`astromol` command for standard output generation. The first-pass CLI lists
registry names, generates individual standard figures, table groups, and slide
decks by stable name, and delegates complete bundle generation through
`astromol outputs`.

Current commands:

- `astromol list`
- `astromol figure cumulative_detections --view current --output cumulative_detections.pdf`
- `astromol table ism_tables --view 2026 --output-dir build/tables`
- `astromol slide ism_molecule_slide --view current --output-dir build/slides`
- `astromol outputs --output-dir build/astromol_outputs --view current --formats png pdf`

Future CLI polish can add shorter aliases and a validation subcommand, but the
core standard-output generation path is no longer blocked.

### 9. Small Documentation Bugs Were Fixed

The README quick-load command now points at the real test path:

```bash
python -m pytest tests/test_load.py
```

It also includes a direct `Database()` load snippet for users who want a quick
interactive check.

### 10. Matplotlib Deprecation Warnings Were Cleaned Up

The horizontal boxplot helper now uses `orientation="horizontal"` instead of
the deprecated `vert=False` argument. This removes Matplotlib deprecation noise
from the production box/strip figures.

## Suggested Order Of Attack

### Current Curation Phase

1. Resolve the remaining inherited `*` dipole-moment placeholders, preserving
   explicit references and notes for any unresolved/null replacements.
2. Stage new molecule and detection updates through YAML templates in
   `curation/staging/`. Codex may draft YAML from maintainer-provided data, but
   the maintainer should review the staged record before `--apply`.
3. Refresh the committed production data baseline after each accepted curation
   batch with `python scripts/update_data_baseline.py`, then inspect
   `tests/baselines/production_data.json` so count and warning changes are
   explicit.
4. Run focused validation after each applied curation batch:
   `python -m astromol.validation` and
   `python -m pytest tests/test_load.py tests/test_validation.py`.
5. Keep broad code changes out of the curation path unless the existing
   workflow blocks a real data update.

### Immediate Cleanup

No immediate cleanup items remain from this review pass.

### Before Heavy 2026 Data Expansion

1. Fix source/telescope template drift and add template/schema consistency
   tests. **Done:** `curation/templates/source.yaml` and
   `curation/templates/telescope.yaml` now expose history consistently, and
   `tests/test_curation_templates.py`
   now checks that each curator template exposes the same production field set
   known to `scripts/stage_records.py`, while allowing staging-only `_...`
   metadata fields.
2. Add duplicate source/telescope nick checks in `Database`. **Done:**
   `Database()` now rejects duplicate source and telescope nicks, and
   `tests/test_database.py` covers both cases.
3. Harden `DipoleMoment.total` around nonnumeric placeholders. **Done:**
   nonnumeric populated components now make `.total` return `None`, with
   coverage in `tests/test_models.py`.
4. Replace `boxplot(vert=False)` with `orientation="horizontal"`.
   **Done:** the shared horizontal boxplot helper now uses the non-deprecated
   Matplotlib argument.
5. Add focused unit tests for model properties and census-view boundaries.
   **Done:** `tests/test_models.py` now covers isotope mass fallback,
   charge/radical inference, DU limits, and kappa behavior;
   `tests/test_census_view_boundaries.py` covers accepted-vs-introduced
   census selection, tentative/disputed inclusion, isotopologue filtering, and
   context molecule selection.
3. Add staging-script tests using small temporary YAML inputs. **Done:**
   `tests/test_stage_records.py` exercises preview-only staging, `--apply`,
   default/history normalization, staging-only metadata reporting, and bad
   reference rejection without touching production JSON.
4. Decide whether schema definitions should become generated templates.
   **Decided/deferred:** do not pursue generated templates during the current
   curation phase. Keep hand-written curator templates protected by
   consistency tests, and revisit only if repeated data-entry work reveals a
   concrete maintenance problem.

### Before Public PyPI Release

1. Resolve PyPI ownership/name access. **Done:** project ownership for the
   `astromol` PyPI package has been transferred to Brett.
2. Add license metadata and a repository license file. **Done:** the project
   now declares the MIT License in `pyproject.toml` and includes a repository
   `LICENSE` file.
3. Add CI package build/install verification. **Done:** `.github/workflows/package.yml`
   builds the sdist and wheel, runs `twine check`, installs the wheel into a
   clean environment outside the source tree, loads `Database()`, runs
   `astromol-validate`, and smoke-tests `astromol-generate-outputs`.
4. Update Colab install cells from `refactor` branch installs to PyPI installs.
   **Deferred until release:** keep branch-based installs while the refactor
   branch remains the live test target; update to PyPI installs immediately
   after the new package release is published.
5. Confirm README/Read the Docs no longer describe the branch as pre-release
   unless that is still intentional. **Deferred until release/default-branch
   transition:** current pre-release language is accurate until the refactor is
   published and promoted.

### After API Stabilization

1. Split `astromol.figures` into focused implementation modules with stable
   re-exports. **Started:** `astromol.figures` is now a package that preserves
   the public import path through `astromol.figures.__init__`, with the existing
   implementation isolated in `astromol.figures._core` and shared visual style
   constants moved to `astromol.figures.style`. Future passes can move figure
   families out of `_core` incrementally without changing user imports.
2. Add an output registry that powers docs, notebooks, and future CLI commands.
   **Done:** `astromol.registry` now exposes canonical figure, table, and slide
   specs with stable names, labels, descriptions, and generation callables.
   `astromol.outputs` consumes the registry for standard GitHub Pages bundles,
   and docs/notebooks now point users/developers to the registry as the
   standard-output inventory.
3. Introduce a first-pass public CLI for standard figures, tables, and slides.
   **Done:** `astromol.cli` installs an `astromol` command with registry-backed
   `list`, `figure`, `table`, `slide`, and `outputs` subcommands. Focused tests
   in `tests/test_cli.py` cover selective generation with fake registry specs.
4. Polish the generated-output GitHub Pages landing page. The current page is
   functional and exposes the right products, but it should eventually get a
   cleaner layout, clearer visual hierarchy, and more polished project branding.
   **Done:** the generated `index.html` now presents a structured landing page
   with a hero summary, featured primary downloads, grouped registry-driven
   inventory sections, and explicit file-path labels for alternate formats and
   manuscript fragments.

## Verification During This Review

- Checkpoint commit `6783be5` was pushed before the review began.
- `python -m astromol.validation` passes with 0 errors and the 3 known
  dipole-placeholder warnings.
- Figure package split verification:
  - `from astromol.figures import ...` smoke test passed.
  - `python -m pytest tests/test_regression_scripts.py` passed
    (29 tests in 143.94 s).
  - `python -m pytest tests/test_load.py tests/test_validation.py` passed.
  - `python -m pytest` passed (52 tests in 152.99 s).
  - `python -m sphinx -W -b html docs /private/tmp/astromol_docs_figures_pkg_check`
    passed.
  - Clean `python -m build --sdist --wheel` passed, and the wheel contains the
    new `astromol/figures/` package without the stale `astromol/figures.py`
    module.
- Output registry verification:
  - registry import smoke test reports 22 figure specs, 9 table specs, and
    2 slide specs.
  - `python -m pytest tests/test_registry.py` passed.
  - `tests/test_outputs.py` now adds a direct `generate_standard_outputs()`
    integration check with fake registry specs, covering registry-driven figure,
    table, slide, report, zip, and `index.html` generation without depending on
    the full scientific rendering stack.
  - the generated-output landing page now exposes featured downloads plus
    grouped bundle/slide/PNG/PDF/table/report sections, with the new structure
    covered by `tests/test_outputs.py`.
  - `python -m astromol.outputs --output-dir /private/tmp/astromol_registry_full_outputs_check --view 2026 --formats png`
    passed and generated 55 products through the registry.
  - `python -m pytest` passed after the registry change (56 tests in 166.58 s).
  - `python -m sphinx -W -b html docs /private/tmp/astromol_docs_registry_check`
    passed.
  - Clean `python -m build --sdist --wheel` passed; the wheel includes
    `astromol/registry.py` and the new `astromol/figures/` package, with no
    stale `astromol/figures.py` module.
- Working tree was clean before this review artifact was added.
