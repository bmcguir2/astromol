# astromol Codebase Review

Date: 2026-05-13
Baseline commit: `6783be5` (`Polish package release metadata`)

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

The main remaining risk is not that the package is conceptually wrong. It is
that several working systems have grown quickly and now encode the same schema,
visual style, or output list in more than one place. That creates drift risk as
the database expands for the 2026 census and as the public API becomes more
stable.

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
  reports 0 errors and 5 known dipole-placeholder warnings; the pytest
  regression harness exercises table, figure, and slide outputs.

## Priority Findings

### 1. Curation Schema Drift Is The Biggest Practical Risk

The curation schema is currently repeated in several places:

- dataclasses and allowed values in `astromol/models.py`
- staging field order/defaults/required fields in `scripts/stage_records.py`
- templates in `curation/templates/*.yaml`
- prose schema in `SPEC.md` and docs
- validation in `astromol/validation.py`

There is already mild drift:

- `scripts/stage_records.py` includes `history` for source and telescope
  records, but `curation/templates/source.yaml` and
  `curation/templates/telescope.yaml` do not expose a history block.
- `curation/templates/telescope.yaml` refers to `TELESCOPE_TYPES`, but no such
  constant exists in `astromol.models`.

Recommended direction: introduce a single schema source for curator-facing
fields, defaults, required markers, and allowed values. It does not need to be
heavyweight. A small `astromol/schema.py` or `curation/schema.py` that generates
templates and feeds `stage_records.py` would remove most of this drift.

Short-term corrective action:

- add source/telescope `history` blocks to the templates or deliberately remove
  them from the exposed template policy;
- fix the `TELESCOPE_TYPES` wording;
- add a small test that checks template keys against `stage_records.py`
  `FIELDS`.

### 2. Source And Telescope Duplicate Keys Can Be Silently Overwritten

`Database._load_molecules` and `_load_detections` reject duplicate primary
keys, but `_load_sources` and `_load_telescopes` currently assign directly into
dictionaries without a duplicate check. A duplicate `nick` in JSON would keep
only the later record in the loaded database.

Recommended action: add duplicate checks for source and telescope nicks during
load, matching the molecule/detection behavior, and add validation tests.

### 3. Known Dipole Placeholders Are Warnings But Can Still Break API Calls

The validator correctly reports the five inherited `*` dipole placeholders as
warnings, but `DipoleMoment.total` attempts numeric exponentiation and raises a
`TypeError` if a user accesses it for those records.

Current warning records:

- `mol:AlCl`
- `mol:CP`
- `mol:SO+`
- `mol:MgCN`
- `mol:HNCS`

Recommended action: keep the manuscript-note reminder to resolve the values,
but also harden `DipoleMoment.total` so nonnumeric placeholders return `None`
or raise a clearer domain-specific error. Since these are currently allowed as
known warnings, the public API should not crash with a raw `TypeError`.

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

### 8. Public CLI Is Minimal

The package currently exposes `astromol-validate`. That is enough for the
current refactor, especially with Colab notebooks, but public users will likely
expect command-line entry points for common outputs.

Recommended future CLI:

- `astromol validate`
- `astromol figure cumulative --view current --format pdf`
- `astromol slide ism --view current`
- `astromol slide ppd --view current`
- `astromol table ism --view 2026`

This should wait until the function registry and figure/table/slide naming are
stable.

### 9. Small Documentation Bugs Remain

The README quick-load command says:

```bash
python test_load.py
```

There is no root-level `test_load.py`; the test lives at `tests/test_load.py`.
Replace this with either:

```bash
python -m pytest tests/test_load.py
```

or a direct `python -c` load snippet. This is low risk but user-facing.

### 10. Matplotlib Deprecation Warnings Are Known

The pytest suite currently passes with Matplotlib deprecation warnings from
`boxplot(vert=False)`. Updating those calls to `orientation="horizontal"` will
remove noise before release.

## Suggested Order Of Attack

### Immediate Cleanup

1. Fix the README quick-load command.
2. Fix source/telescope template drift and the `TELESCOPE_TYPES` comment.
3. Add duplicate source/telescope nick checks in `Database`.
4. Harden `DipoleMoment.total` around nonnumeric placeholders.
5. Replace `boxplot(vert=False)` with `orientation="horizontal"`.

### Before Heavy 2026 Data Expansion

1. Add template/schema consistency tests. **Done:** `tests/test_curation_templates.py`
   now checks that each curator template exposes the same production field set
   known to `scripts/stage_records.py`, while allowing staging-only `_...`
   metadata fields.
2. Add focused unit tests for model properties and census-view boundaries.
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
   **Decided/deferred:** keep hand-written curator templates for now, protected
   by consistency tests. Revisit generated templates after the 2026 data model
   is less fluid, so the generator does not harden a schema that is still
   actively changing.

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
   re-exports.
2. Add an output registry that powers docs, notebooks, and future CLI commands.
3. Introduce a first-pass public CLI for standard figures, tables, and slides.
4. Polish the generated-output GitHub Pages landing page. The current page is
   functional and exposes the right products, but it should eventually get a
   cleaner layout, clearer visual hierarchy, and more polished project branding.

## Verification During This Review

- Checkpoint commit `6783be5` was pushed before the review began.
- `python -m astromol.validation` passes with 0 errors and the 5 known
  dipole-placeholder warnings.
- Working tree was clean before this review artifact was added.
