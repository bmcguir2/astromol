# Curation Workflow

New records should be staged through YAML before they are applied to production
JSON.

The current curation workflow is a source workflow (Git checkout or extracted source archive) rather than an
installed package command. From the repository root, install the contributor
dependencies first:

```bash
python -m pip install -e ".[dev]"
```

1. Copy one or more templates from `curation/templates/`.
2. Fill in a staging file under `curation/staging/`.
3. Generate a preview and review the report.
4. Apply the staged records after human approval.
5. Remove the staging YAML after the production update is committed.

Generate a preview:

```bash
python scripts/stage_records.py --staging curation/staging/example.yaml
```

Apply approved staged records:

```bash
python scripts/stage_records.py --staging curation/staging/example.yaml --apply
```

The staging templates intentionally expose the full available field set,
including optional fields, so curators do not need to remember schema details
from memory.

After an approved production data change, refresh the committed data baseline
and inspect the diff:

```bash
python scripts/update_data_baseline.py
python scripts/check_curation.py
```

`tests/baselines/production_data.json` stores expected production inventory
counts, known validation warnings, and curation-sensitive current regression
expectations. Regenerate it only after the data change has been reviewed and
accepted. Add new scientific membership/count expectations to the generated baseline;
keep evolving layout and formatting expectations in focused invariant tests.

`scripts/check_curation.py` runs production-data validation plus the load,
validation-baseline, and output regression tests that are sensitive to new
secure molecules and detections.

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

Scientific data are human curated and human verified. LLM assistance may be used
for code, staging support, or mechanical transformations, but production data
are not accepted without curator review.

## Reviewed Preview And Recovery

Apply requires a successful preview for the same batch name. Review the report,
full preview JSON, and proposed baseline before passing `--apply`. Any change to
YAML, production files, bibliography, preview artifacts, curation code, or relevant
dependencies requires a fresh preview and review. Applying later preserves the
preview date. All proposed files and the baseline are prepared before production
replacement; ordinary failures restore originals. If execution was interrupted,
run `python scripts/stage_records.py --recover <name>` before retrying. Recovery
refuses to overwrite files edited after the interrupted apply.

New working records receive dated history without assuming a publication-year
census label. `current` includes accepted records and optional tracked tentative
or disputed records; `2026` is a compatibility alias. Published 2018/2021
membership continues to use explicit history census labels.

The baseline records reviewed membership, counts, numerical summaries, and known
warnings. Its analysis endpoint is the latest secure detection year, keeping
snapshots independent of the day tests run. Current rendering tests check unique
membership, complete grouping, status markers, valid table structure, and overflow
instead of recording fonts or formatted labels in this baseline. A generated
snapshot is a change-review aid; independent fixtures also test selection rules.

Commit-mode cleanup refuses unrelated staged paths before deleting anything and
restores review files if the commit fails. A push failure after a successful
commit leaves the committed cleanup intact so the push can be retried.
If cleanup itself was interrupted, restore its saved review artifacts with
`python scripts/cleanup_stage.py --name <name> --recover` before retrying.
