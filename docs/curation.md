# Curation Workflow

Stage additions and updates in YAML so the maintainer can review the complete
proposed records before they enter production.

The curation scripts run from a Git checkout or extracted source archive.
From its root directory, install the contributor tools:

```bash
python -m pip install -e ".[dev]"
```

1. For new records, copy templates from `curation/templates/`. For existing
   records, generate update templates with `--prepare-update` (see below).
2. Fill in a staging file under `curation/staging/`.
3. Generate a preview and review the report.
4. Apply the staged records after human approval.
5. Inspect the generated baseline diff and run the curation check.
6. Use the cleanup helper when the batch is finished.

Generate a preview:

```bash
python scripts/stage_records.py --staging curation/staging/example.yaml
```

Apply approved staged records:

```bash
python scripts/stage_records.py --staging curation/staging/example.yaml --apply
```

The templates show every field, including optional fields, so you can fill in
the record without remembering the schema. Leave unknown optional values blank.
For an existing record, generate a template containing its current values:

```bash
python scripts/stage_records.py \
  --prepare-update molecule mol:EXAMPLE \
  --output curation/staging/example_update.yaml
```

Add an `_update_summary` and edit the fields that need to change. The template's
`_base_digest` lets staging check whether the production record has changed
since you prepared the update.

Apply refreshes `tests/baselines/production_data.json` automatically. Inspect
that diff, then run:

```bash
python scripts/check_curation.py
```

The baseline records expected membership, counts, numerical summaries, and
accepted validation warnings. New expectations for scientific membership or
counts belong there; layout and formatting checks belong in focused tests.
Use `python scripts/update_data_baseline.py` for a manual baseline refresh
outside the staging workflow, after the underlying change has been reviewed.

`scripts/check_curation.py` validates the data and runs the load, baseline, and
output regression tests. New records can change counts and the space needed
by figures, tables, and slides, so validation alone is not enough.

When the batch is finished, remove its review files with:

```bash
python scripts/cleanup_stage.py --name example
```

References are maintained in Zotero and exported to
`astromol/data/references.bib`. Do not hand-edit the BibTeX file; export the
updated collection and commit the resulting file with the relevant data change.

## Project-Computed Values

Some curated values may be calculated as part of this work when no suitable
literature value is available. Keep the supporting notebook under
`docs/calculations/` and commit it with the data update. Do not commit local
`.ipynb_checkpoints/` directories.

Structured molecule reference fields resolve to BibTeX citekeys, so use those
fields for the software, method, basis-set, or literature references that
support the calculation. Record the local notebook path in the relevant note or
history entry. If the calculation is later described in a citable publication,
add that citekey through the normal Zotero export workflow.

The maintainer verifies all scientific data before apply. AI assistance can
help with code, staging, and transformations of supplied records.

## Reviewed Preview And Recovery

Apply uses the reviewed preview for that batch and preserves its date. Changes
to the YAML, production data, bibliography, review files, curation code, or
relevant dependencies require a fresh preview and review.

The script prepares all proposed files and the baseline before replacing
production data. Ordinary failures restore the originals. After an interrupted
apply, run `python scripts/stage_records.py --recover <name>` before retrying.
Recovery refuses to overwrite files edited after the interruption.

New working records receive dated history without an assumed census year.
`current` selects accepted records; tentative and disputed records can be
requested explicitly. `2026` is an alias for `current`. Published 2018/2021
membership still uses the census labels in record history.

The baseline uses the latest secure detection year as its analysis endpoint,
so its expected values do not change with the day tests run. Rendering tests
check membership, grouping, status markers, table structure, and overflow.
Separate fixtures test selection rules independently of the generated baseline.

When asked to commit, cleanup checks for unrelated staged files before deleting
anything and restores the review files if the commit fails. If the commit
succeeds but the push fails, the committed cleanup stays in place so you can
retry the push.
If cleanup itself was interrupted, restore its saved review artifacts with
`python scripts/cleanup_stage.py --name <name> --recover` before retrying.
