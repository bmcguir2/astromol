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
2026 regression expectations are also checked against the committed baseline at
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
