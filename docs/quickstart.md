# Quickstart

Load the production database:

```python
from astromol.database import Database

db = Database()
print(len(db.molecules), len(db.detections))
```

Create a frozen census view:

```python
from astromol.census import CensusView

view_2021 = CensusView.for_census(db, "2021")
view_2026 = CensusView.for_census(db, "2026")
current = CensusView.current(db)
```

Query molecules in a context:

```python
ism_molecules = view_2026.ism_molecules()
ppd_molecules_with_isotopologues = view_2026.ppd_molecules(
    include_isotopologues=True
)
```

By default, context views exclude isotopologues and include only secure
accepted detections. Tentative and disputed detections can be requested
explicitly:

```python
exgal_with_tentatives = view_2026.exgal_detections(include_tentative=True)
```

Run production-data validation:

```bash
astromol-validate
```

Generate a molecule slide:

```python
from pathlib import Path

from astromol.slides import write_molecule_slide

write_molecule_slide(view_2026, Path("astro_molecules.pptx"), profile="balanced")
```
