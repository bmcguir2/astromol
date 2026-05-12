# astromol documentation

`astromol` is a Python package and curated database for astronomical molecule
detections. It stores molecule, detection, source, telescope, and reference
metadata in version-controlled files, then exposes those records through a
small Python API for validation, curation, census views, manuscript tables,
figures, and PowerPoint slides.

This documentation describes the active refactor branch. The data model,
curation workflow, and output-generation tools are usable, but the public API is
still pre-release.

```{toctree}
:maxdepth: 2
:caption: User Guide

installation
quickstart
faq
data-model
curation
generation
census-outputs
notebooks
validation
```

```{toctree}
:maxdepth: 2
:caption: Reference

api
```

## Project Records

The repository also contains durable project records that are useful while the
refactor is still active:

- `SPEC.md`: current architecture and schema specification
- `GENERATION_MIGRATION.md`: table, figure, and slide migration audit trail
- `MANUSCRIPT_NOTES_2026.md`: writing-time reminders for the 2026 census paper
