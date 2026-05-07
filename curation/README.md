# astromol Curation Workflow

This directory contains curator-facing templates for staging new database
records without hand-editing the production JSON files directly.

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

After reviewing the preview report, apply with:

```bash
python scripts/stage_records.py --staging curation/staging/example.yaml --apply
```

## Record Kinds

Each staged record has a `kind`:

- `molecule`
- `detection`
- `source`
- `telescope`

Fields beginning with `_` are staging-only notes. They are allowed in YAML and
reported when useful, but they are not written into production JSON.

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
`history.accepted` block. For a molecule that should be tracked but is not yet
accepted as confirmed, set `history.accepted: null` and set
`history.introduced.context` to `tentative` or `disputed`.

Detection records require a stable `id`. Use the default format
`det:<molecule-label-without-mol-prefix>:<context>:<year>`, for example
`det:CH3CH2CCH:ism-csm:2021`. Add a final qualifier only when that default
would collide with another detection. Keep mutable status words such as
`tentative` out of IDs.

## YAML Value Notes

YAML will sometimes reinterpret unquoted values. Quote values that must remain
strings when they contain colon-delimited coordinates, leading signs, or other
special syntax.

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
