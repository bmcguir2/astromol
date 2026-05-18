# Curation Workflow

New records should be staged through YAML before they are applied to production
JSON.

The current curation workflow is a source-checkout workflow rather than an
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
