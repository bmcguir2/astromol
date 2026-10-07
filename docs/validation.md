# Validation And Tests

Run semantic production-data validation before committing data updates:

```bash
astromol-validate
```

or from a source checkout:

```bash
python -m astromol.validation
```

Validation errors should block release and data-update commits. Validation
warnings are allowed for known unresolved curation follow-ups, such as inherited
legacy `*` dipole placeholders.

Production inventory counts, known validation warnings, and curation-sensitive
current regression expectations are also checked against the committed baseline at
`tests/baselines/production_data.json`. Refresh that file only after a reviewed
data update:

```bash
python scripts/update_data_baseline.py
python scripts/check_curation.py
```

Do not repair curation-driven assertion failures by editing scattered numeric
literals in regression scripts. Add the generated expectation to
`scripts/update_data_baseline.py`, regenerate
`tests/baselines/production_data.json`, and inspect the baseline diff with the
curation change.

`scripts/check_curation.py` expands to production-data validation plus the
load, validation-baseline, and output regression tests. Use it after applying
new molecule or detection records.

Run the regression suite with:

```bash
python -m pytest
```

The current pytest suite includes:

- native pytest checks for database loading and validation
- a wrapped regression-script harness for the migration-era table, figure, and
  slide checks

The wrapper is intentional for now: it preserves the audited migration checks
while making them visible to pytest. Once the public API, documentation
examples, and CI workflow stabilize, the highest-value checks should be
converted incrementally into conventional unit and integration tests.

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
