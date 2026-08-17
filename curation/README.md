# astromol Curation Workflow

This directory contains curator-facing templates for staging database
additions and updates without hand-editing production JSON directly.

## Basic Workflow

1. Copy a template from `curation/templates/`.
2. Paste it into a YAML file under `curation/staging/`.
3. Fill in the fields you know. Fields marked `# REQUIRED` must be filled in.
   Optional fields can be left blank; the staging script prunes empty template
   sections and fills production defaults.
4. Run:

```bash
python scripts/stage_records.py --staging curation/staging/example.yaml
```

This writes preview JSON and a report under `astromol/data/`. It does not
modify production data unless `--apply` is passed.

`stage_records.py` only reads the file paths passed with `--staging`; it does
not automatically process every YAML file in `curation/staging/`. Use
`curation/staging/` for active records you are currently testing, and
`curation/todo/` for backlog notes or generated to-do lists.

The preview report includes a `Generated Count Updates` section showing any
production inventory or regression-count baseline values that would change if
the staged records are applied. After reviewing the preview report, apply with:

```bash
python scripts/stage_records.py --staging curation/staging/example.yaml --apply
```

Applying staged records refreshes
`tests/baselines/production_data.json` automatically. Do not manually edit
generated count expectations in regression scripts; update the curated data and
let the staging workflow regenerate the baseline.

### Updating Existing Records

Each staging record has an `operation` field. `operation: add` is the default
and remains backward compatible with older staging YAML. Generate a full-field
template before using `operation: update`:

```bash
python scripts/stage_records.py \
  --prepare-update molecule mol:EXAMPLE \
  --output curation/staging/example_molecule_update.yaml
```

The generated block includes the current production record and a staging-only
`_base_digest`. Fill in `_update_summary`, make the intended changes, and leave
the other full-template fields visible and unchanged. Staging fails if the
production record has changed since the template was generated. The workflow
derives changed field paths, refreshes `history.last_modified`, and appends an
update event. Use `_event_kind: corrected` for corrections.

Additions and updates may share one batch. Detection relationship reciprocals
are mechanical: declaring `confirms`, `disputes`, or `supersedes` derives and
reports the target record's reciprocal update in the merged preview. Scientific
changes such as molecule promotion remain explicit `operation: update`
records. The complete merged preview must load and pass semantic validation
before `--apply` writes any production file.

Successful staging runs also write a manifest under `astromol/data/`, for
example `astromol/data/example_stage_manifest.json`. Use that manifest with the
cleanup helper after the curation batch is finished:

```bash
python scripts/cleanup_stage.py --name example
```

This removes the staged YAML input, preview JSON artifacts, stage report, and
the manifest itself. When committing, the helper also auto-stages the most
common curation sidecar files if they are modified:

- `astromol/data/references.bib`
- `tests/baselines/production_data.json`

Applied manifests contain hashes for every touched production file. Cleanup
rejects post-apply drift. When `--commit-message` is used, cleanup also runs
`python scripts/check_curation.py` before deleting review artifacts or
committing. `--skip-verification` is an explicit recovery-only escape hatch;
it does not bypass manifest hash verification.

For additional files beyond those defaults, use `--include`. To create a
commit and push it:

```bash
python scripts/cleanup_stage.py \
  --name example \
  --commit-message "Add Example curation batch" \
  --push \
  --close-issue 123
```

`--push` is intentionally separate from `--commit-message` so remote updates
remain an explicit choice. `--close-issue` may be passed more than once and
runs only after a successful push; it uses the GitHub CLI to close each issue
with a comment linking to the commit.

After any approved production data change, review the generated baseline diff
and run the local curation verification check:

```bash
python scripts/check_curation.py
```

`tests/baselines/production_data.json` records the expected production
inventory counts, curation-sensitive regression counts, and known validation
warnings. It should change only when the underlying curated data or accepted
warning backlog changes. The curation check runs production-data validation
plus the load, validation-baseline, and output regression tests that catch
2026 table, figure, and slide expectation drift. The
`python scripts/update_data_baseline.py` command remains available for manual
baseline refreshes outside the staging workflow.

## Record Kinds

Each staged record has a `kind`:

- `molecule`
- `detection`
- `source`
- `telescope`

Fields beginning with `_` are staging-only notes. They are allowed in YAML and
reported when useful, but they are not written into production JSON.

## Project-Computed Values

When a database value is calculated as part of this project, keep the
supporting notebook under `docs/calculations/` and commit it with the data
change. The structured `refs` fields should continue to contain BibTeX citekeys
for software, methods, basis sets, or literature sources that support the
calculation. Record the local notebook path in an appropriate note or history
summary. Do not track local `.ipynb_checkpoints/` directories.

## Generated Staging Files

When Codex or another helper generates staging YAML for curator review, it
should include the full field set from the relevant template for every staged
record. Do not emit compact records containing only fields that were obvious
from the source text. Curators need the unused optional fields visible so they
can fill in laboratory constants, dipole moments, identifiers, notes,
relationships, and other metadata during review.

Leave unknown optional values blank. The staging script prunes blank template
values and fills production defaults when previewing or applying records.

History metadata is generated automatically for staged records. By default the
script adds `history.introduced.date`, `history.last_modified`, and an initial
dated `added` event using the staging run date. New molecule records also
default to `history.introduced.context: confirmed` and a current-census
`history.accepted` block. Secure detection records likewise default to a
current-census `history.accepted` block for their detection context. Tentative
or disputed detection records default to `history.accepted: null`. For any
record that should be tracked but is not yet accepted as confirmed, set
`history.accepted: null` and set `history.introduced.context` to `tentative` or
`disputed`.

Detection records require a stable `id`. Use the default format
`det:<molecule-label-without-mol-prefix>:<context>:<year>`, for example
`det:CH3CH2CCH:ism-csm:2021`. Add a final qualifier only when that default
would collide with another detection. Keep mutable status words such as
`tentative` out of IDs.

## YAML Value Notes

YAML will sometimes reinterpret unquoted values. Quote values that must remain
strings when they contain colon-delimited coordinates, leading signs, or other
special syntax.

Quote ISO date strings in staging YAML:

```yaml
date: "2025-01-07"
last_modified: "2025-01-07"
```

Do not leave bare `YYYY-MM-DD` values unquoted in fields such as
`history.introduced.date`, `history.accepted.date`, or `history.last_modified`.
YAML may parse them as native dates, and `stage_records.py` writes JSON
previews that expect ordinary strings.

Use quoted strings for source coordinates:

```yaml
ra: "18:53:18.5"
dec: "+01:14:59"
```

Do not use sexagesimal strings for telescope site coordinates. Telescope
`latitude` and `longitude` are decimal degrees and should be entered as numbers:

```yaml
latitude: 38.433056
longitude: -79.839722
```

## Schema Maintenance

When `Molecule`, `Detection`, `Source`, or `Telescope` changes in
`astromol/models.py`, update all of the following in the same change:

- the relevant template in `curation/templates/`
- the defaults and field order in `scripts/stage_records.py`
- the relevant data-model section in `SPEC.md`

This keeps the curator-facing workflow synchronized with the production schema.
