# Output Generation

Each output separates selecting and calculating the data from rendering it.
That lets tests check the records and numerical results before a figure, table,
or slide is drawn.

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
brand colors only when they serve a clear visual purpose. Figure writers save
PNG output, and any rasterized artists embedded in PDF output, at 300 DPI.

`astromol.registry` lists the standard outputs used by the CLI and bundle
generator. To see what is available, or to add a standard output, start there:

```python
from astromol.registry import FIGURE_OUTPUTS, SLIDE_OUTPUTS, TABLE_OUTPUTS

for spec in FIGURE_OUTPUTS:
    print(spec.name, "-", spec.label)
```

`astromol-generate-outputs` also writes a static `index.html` landing page for
GitHub Pages. It links the complete bundle, cumulative detections figure, and
two slide decks first, followed by the other figures. The bundle also contains
LaTeX fragments and slide layout reports.

The public `astromol` command can list and generate individual registry
products:

```bash
astromol list all
astromol figure cumulative_detections --view current --output cumulative_detections.pdf
astromol table ism_tables --view current --output-dir build/tables
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
view = CensusView.current(db)

write_molecule_slide(view, Path("astro_molecules_current.pptx"), profile="balanced")
write_ppd_detection_slide(view, Path("ppd_molecules_current.pptx"))
```

For examples you can copy and run for each census output, see
[](census-outputs.md).

## Reproducible Bundles

Use a new or empty directory for the first bundle. The generator validates
inputs and builds in a temporary directory before replacing an existing bundle.
It will replace a nonempty directory only if it is marked as an astromol bundle.
A failed build preserves the previous bundle. `--no-clean`
retains other files in an owned directory, while the zip includes only products
from this run, the landing page, and `manifest.json`.

`manifest.json` records the source revision and uncommitted changes when Git
metadata are available, along with code/data hashes, dependency versions,
settings, output hashes, and the analysis end year. An installed wheel without
Git metadata records a null revision but still records code hashes. Current
analyses default to the calendar year; choose an end year explicitly when you
need to repeat an analysis:

```bash
astromol outputs --view current --end-year 2026 --output-dir build/paper_outputs
```

This endpoint controls rate calculations and plot ranges, not census membership.
Slide update dates include selected molecule, detection, source, and telescope
history. Checkout slides show the Git revision even with installed development
package metadata.

During the refactor, only `refactor` pushes publish Pages, after data validation,
tests, documentation, and output generation succeed. Other branches produce
review artifacts. Read the Docs and Colab continue following `refactor`.
