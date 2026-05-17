# Output Generation

Output-generation helpers are organized around reusable data builders and writer
functions. The data builders are tested separately from visual rendering so
scientific selections can be verified without relying only on visual inspection.

## Tables And LaTeX

LaTeX/table helpers live in `astromol.latex`. They consume `CensusView`
instances and generate manuscript-ready fragments for historical reproduction
or current-census output.

## Figures

Figure helpers live in `astromol.figures`. Migrated figures generally provide:

- a data-builder function
- a plotting function
- a writer function

The current production style uses readable manuscript-scale typography,
color-blind friendlier palettes when multiple categories are shown, and project
brand colors only when they serve a clear visual purpose.

The standard output inventory lives in `astromol.registry`. The registry is the
canonical list used by the generated-output bundle and should be the first place
to look when adding a new production figure, table, or slide:

```python
from astromol.registry import FIGURE_OUTPUTS, SLIDE_OUTPUTS, TABLE_OUTPUTS

for spec in FIGURE_OUTPUTS:
    print(spec.name, "-", spec.label)
```

`astromol-generate-outputs` also writes a static `index.html` landing page for
GitHub Pages. The page surfaces the standard bundle and key slide/figure
downloads first, then presents public-facing figure download cards. LaTeX
fragments and layout diagnostics remain available in the complete bundle.

The public `astromol` command can list and generate individual registry
products:

```bash
astromol list all
astromol figure cumulative_detections --view current --output cumulative_detections.pdf
astromol table ism_tables --view 2026 --output-dir build/tables
astromol slide ppd_detection_slide --view current --output-dir build/slides
```

## Slides

PowerPoint helpers live in `astromol.slides`.

The ISM/CSM molecule slide supports:

- `profile="legacy"` for 2021 layout reproduction
- `profile="balanced"` for the dynamic current-census layout

The PPD detection slide includes isotopologues by default and uses a compact
layout for the smaller PPD inventory.

Example:

```python
from pathlib import Path

from astromol.census import CensusView
from astromol.database import Database
from astromol.slides import write_molecule_slide, write_ppd_detection_slide

db = Database()
view = CensusView.for_census(db, "2026")

write_molecule_slide(view, Path("astro_molecules_2026.pptx"), profile="balanced")
write_ppd_detection_slide(view, Path("ppd_molecules_2026.pptx"))
```

For detailed, copy-pasteable recipes covering every migrated census table,
figure, and slide product, see [](census-outputs.md).
