# Validation And Tests

Validate production data before committing an update:

```bash
astromol-validate
```

or from a source checkout:

```bash
python -m astromol.validation
```

Resolve validation errors before a release or data-update commit. Warnings
track unresolved curation follow-ups that the maintainer has accepted; the
current production baseline has no accepted warnings.

Tests also compare membership, counts, numerical summaries, and warnings with
`tests/baselines/production_data.json`. Staging refreshes that baseline on apply.
Inspect the diff, then run the curation check:

```bash
python scripts/check_curation.py
```

If a new scientific count needs regression coverage, add its calculation to
`scripts/update_data_baseline.py` so the expectation is generated from the data.
Use that script for manual baseline refreshes outside staging, after review.
Do not patch numeric expectations by hand to make a failed test pass.

The curation check validates the data and runs the load, baseline, and output
regression tests. Use it after applying records: a valid record can still
change a count or make a table or slide overflow.

Run the regression suite with:

```bash
python -m pytest
```

The current pytest suite includes:

- native pytest checks for database loading and validation
- a wrapped regression-script harness for the migration-era table, figure, and
  slide checks

The wrapper lets pytest run the audited migration checks. As the API, examples,
and CI workflow stabilize, those checks can move into focused unit and
integration tests.

## Reviewed Preview And Recovery

Review the full preview and proposed baseline before apply. The script checks
that the reviewed inputs are unchanged and preserves the preview date. Changed
inputs require a fresh preview and review.

Ordinary apply failures restore the originals. After an interruption, run
`python scripts/stage_records.py --recover <name>` before retrying. Recovery
will not overwrite files edited since the interruption. See [](curation.md)
for the full review and cleanup workflow.

The baseline uses the latest secure detection year as its analysis endpoint,
so expected values do not change with the day tests run. Rendering tests check
membership, grouping, status markers, table structure, and overflow. Separate
fixtures check selection rules independently of the generated baseline.
