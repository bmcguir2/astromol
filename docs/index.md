# astromol documentation

`astromol` is a Python package and curated database for astronomical molecule
detections. It stores molecule, detection, source, telescope, and reference
metadata in version-controlled files, then exposes those records through a
small Python API for validation, curation, census views, manuscript tables,
figures, and PowerPoint slides.

This documentation describes the active refactor branch. The data model,
curation workflow, and output-generation tools are usable, but the public API is
still pre-release.

## Latest Generated Products

The latest standard figures and PowerPoint slides are generated automatically
from the current database and published at:

<https://bmcguir2.github.io/astromol/>

Use that page for the latest ISM/CSM cumulative detections figure, the ISM/CSM
molecule slide, the PPD molecule/isotopologue slide, and the complete generated
output bundle. Use the Colab notebooks when you want to customize the view or
regenerate selected outputs interactively.

```{toctree}
:maxdepth: 2
:caption: User Guide

installation
quickstart
faq
data-model
curation
calculations/README
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
