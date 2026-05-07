# astromol

`astromol` is a Python library and curated data set for cataloging astronomical
molecule detections. It stores molecule, detection, source, telescope, and
reference metadata in version-controlled data files, then loads and resolves
those records through a small Python API.

This README describes the active refactor branch. The data model and curation
workflow are usable, but the public package/API documentation is still being
developed before this branch replaces `main`.

## Current Status

- Production data lives in `astromol/data/*.json`.
- Bibliographic metadata lives in `astromol/data/references.bib`.
- `references.bib` is exported from Zotero and should not be hand-edited.
- `SPEC.md` is the source of truth for the current schema and architecture.
- New records should be staged through YAML files in `curation/staging/` before
  being applied to production JSON.

## Quick Load Check

```bash
python test_load.py
```

Or from Python:

```python
from astromol.database import Database

db = Database()
print(len(db.molecules), len(db.detections))
```

`Database()` loads references, telescopes, sources, molecules, and detections;
it also resolves cross-references and validates stable detection IDs.

## Curation Workflow

Copy a template from `curation/templates/` into `curation/staging/`, fill in the
fields you know, then generate a preview:

```bash
python scripts/stage_records.py --staging curation/staging/example.yaml
```

After reviewing the generated report and preview JSON under `astromol/data/`,
apply the staged records:

```bash
python scripts/stage_records.py --staging curation/staging/example.yaml --apply
```

See `curation/README.md` for template notes, YAML quoting rules, and schema
maintenance expectations.

## Reference Workflow

References are maintained in Zotero and exported to
`astromol/data/references.bib`. When new data records use new citekeys, export
the updated Zotero collection and commit the resulting BibTeX file.

Internal citekeys use:

```text
LastName:Year:FirstPage
```

Add lowercase suffixes only when needed to disambiguate duplicate keys.

## Important Files

- `SPEC.md`: architecture and schema specification
- `astromol/models.py`: dataclasses and allowed values
- `astromol/database.py`: JSON/BibTeX loader and cross-reference resolver
- `astromol/data/`: production data files
- `curation/templates/`: curator-facing YAML templates
- `scripts/stage_records.py`: staging, validation, preview, and apply workflow

## Packaging

The long-term goal is to make `astromol` installable from PyPI. Packaging files
and final user-facing API documentation have not yet been added on this branch.
