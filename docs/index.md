# astromol documentation

`astromol` is a Python package and curated database of astronomical molecule
detections. Use it to query the inventory or generate figures, manuscript
tables, and PowerPoint slides. Each output uses the same database and census
selection rules, so you can change the view without rebuilding the analysis
by hand.

These docs describe the active `refactor` branch. The database, curation
workflow, and output tools are usable; the package and public API are still
pre-release.

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

The repository keeps the schema and migration checks in these files:

- `SPEC.md`: current architecture and schema specification
- `GENERATION_MIGRATION.md`: table, figure, and slide migration audit trail
- `MANUSCRIPT_NOTES_2026.md`: writing reminders for the next census paper
