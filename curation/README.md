# astromol Curation Workflow

Use the templates here to stage database additions and updates in YAML.
The preview lets you review the proposed records before changing production
JSON.

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

Additions and updates can share a batch. Declare `confirms`, `disputes`, or
`supersedes` on the scientifically meaningful side; the script adds the target
record's reciprocal link and shows it in the preview. Scientific changes such
as molecule promotion must be explicit `operation: update` records. The
complete preview must load and pass validation before apply.

Successful staging runs also write a manifest under `astromol/data/`, for
example `astromol/data/example_stage_manifest.json`. Use that manifest with the
cleanup helper after the curation batch is finished:

```bash
python scripts/cleanup_stage.py --name example
```

This removes the staging YAML, preview JSON, report, and manifest. When asked
to commit, the helper stages the touched production JSON files and these
supporting files if they have changed:

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

Add `--push` when you also want to push the commit. You can pass `--close-issue`
more than once. After a successful push, the helper uses the GitHub CLI to close
each issue with a comment linking to the commit.

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
current table, figure, and slide expectation drift. The
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

When generating staging YAML for review, include every field from the relevant
template. Keep unused optional fields visible so the curator can fill in
constants, dipole moments, identifiers, notes, and relationships during review.

Leave unknown optional values blank. The staging script prunes blank template
values and fills production defaults when previewing or applying records.

History metadata is generated automatically for staged records. By default the
script adds `history.introduced.date`, `history.last_modified`, and an initial
dated `added` event using the staging run date. New molecule records also
default to `history.introduced.context: confirmed` and a dated
`history.accepted` block, without an assumed census year. Secure detection
records likewise default to a
dated `history.accepted` block for their detection context. Tentative
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

The templates and staging script need to accept the same fields as production.

## Reviewed Preview And Recovery

Review the report, full preview JSON, and proposed baseline before apply.
The script requires the reviewed preview for the same batch and preserves its
date. Changes to the YAML, production data, bibliography, review files, curation
code, or relevant dependencies require a fresh preview and review.

All proposed files and the baseline are prepared before replacing production
data. Ordinary failures restore the originals. After an interrupted apply,
run `python scripts/stage_records.py --recover <name>` before retrying.
Recovery refuses to overwrite files edited after the interruption.

New working records receive dated history without an assumed census year.
`current` selects accepted records; tentative and disputed records can be
requested explicitly. `2026` is an alias for `current`. Published 2018/2021
membership still uses the census labels in record history.

The baseline records reviewed membership, counts, numerical summaries, and
accepted warnings. It uses the latest secure detection year as its analysis
endpoint, so expected values do not change with the day tests run. Rendering
tests check membership, grouping, status markers, table structure, and overflow;
separate fixtures test selection rules independently of the generated baseline.

When asked to commit, cleanup checks for unrelated staged files before deleting
anything and restores the review files if the commit fails. If the commit
succeeds but the push fails, the committed cleanup stays in place so you can
retry the push.
If cleanup itself was interrupted, restore its saved review artifacts with
`python scripts/cleanup_stage.py --name <name> --recover` before retrying.
