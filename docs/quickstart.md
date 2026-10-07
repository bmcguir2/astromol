# Quickstart

Load the production database:

```python
from astromol.database import Database

db = Database()
print(len(db.molecules), len(db.detections))
```

Choose the inventory you want to work with. This creates both a published
2021 membership view and a view of the current database:

```python
from astromol.census import CensusView

view_2021 = CensusView.for_census(db, "2021")
current = CensusView.current(db)
```

Query molecules in an astronomical context, such as ISM/CSM or protoplanetary
disks:

```python
ism_molecules = current.ism_molecules()
ppd_molecules_with_isotopologues = current.ppd_molecules(
    include_isotopologues=True
)
```

Context queries include secure, accepted detections and exclude isotopologues
by default. Request tentative or disputed detections explicitly when you need
them for a comparison or review:

```python
exgal_with_tentatives = current.exgal_detections(include_tentative=True)
```

Run production-data validation:

```bash
astromol-validate
```

Generate a molecule slide:

```python
from pathlib import Path

from astromol.slides import write_molecule_slide

write_molecule_slide(current, Path("astro_molecules.pptx"), profile="balanced")
```
