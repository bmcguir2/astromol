# 2026 Manuscript Notes

This file collects writing-time reminders that arose during the database and
output-generation migration. Keep implementation details in `SPEC.md` and
`GENERATION_MIGRATION.md`; keep manuscript-facing explanations here.

## Terminology And Classification

- Explain the source-type terminology change from `LOS Cloud` to
  `Diffuse Cloud`. The new term is meant to describe diffuse-cloud
  line-of-sight material more directly, while still allowing the source table to
  consolidate these sight lines into a single `Diffuse Cloud` row.

## Data Curation Follow-Ups

- Current work is shifting from large-scale code refactoring to database
  curation for the 2026 census. Prioritize resolving known data warnings, then
  stage new molecule and detection records through the YAML workflow for
  maintainer review before applying them to production JSON.
- Legacy `*` dipole-moment placeholders have been resolved locally for SO+,
  MgCN, and HNCS. Re-check the validation baseline before final dipole-based
  analysis or manuscript statements.
- Project-computed values should be reproducible from tracked notebooks under
  `docs/calculations/`. For example, the CP and SO+ dipole-moment calculations
  are kept in `docs/calculations/cp_dipole_calculations.ipynb` and
  `docs/calculations/sop_dipole_calculations.ipynb`.
- Add a short manuscript section labeled `sec:dipole` that describes the
  project-computed dipole workflow, including PySCF, pyberny, M06-2x, and the
  6-311++G(d,p) basis-set references. Individual molecule notes should point
  to that section rather than repeating the full computational-method sentence.

## Tables And Counts

- The standard ISM/CSM and external-galaxy molecule tables are
  isotopologue-free. Isotopologue-expanded views exist for analysis, but table
  captions or surrounding prose should make clear which convention is being
  used.
- The standard ice table includes tentative detections and marks them with a
  dagger. In particular, OCN- was inserted directly into the 2021 ice table for
  simplicity, but the modern database represents the OCN- ice claim as
  tentative. Explain this status explicitly in the 2026 text.
- Extragalactic tables include tentative detections by default and mark them
  with a dagger. PPD and exoplanet tables exclude tentative/disputed detections
  unless a specific expanded review view is being produced.

## Rate Analysis

- Use true `R^2` values in the detection-rate-by-atoms table. The 2021 paper's
  table header said `R^2`, but the values printed there were Pearson `R`
  values from `scipy.stats.linregress`. The 2026 manuscript should note this
  correction if comparing directly to the previous paper.
- The cumulative-detections figure no longer labels individual facility
  commissioning dates in the 2026/current view. Legacy 2018 and 2021 views keep
  those annotations for reproduction, but the 2026 figure should discuss the
  removal if comparing visually to previous census figures.
- The cumulative-by-atoms line figure now has a stacked companion figure. The
  line view preserves continuity with the 2021 figure, while the stacked view
  is intended to show how atom-count classes contribute to total inventory
  growth. These should be placed side-by-side as two panels in a single figure;
  the generated PDFs use matched figure size and axes geometry for that layout.
- The 2026/current cumulative-by-atoms panels use a solid color-blind-friendly
  palette. Legacy 2018/2021 reproduction views preserve the original palette.
  Dashed lines and hatches were tested but rejected because they reduced
  readability in the dense line and stacked panels.
- Explicitly mention in the 2026 manuscript that the figure palettes were
  updated to be more color-blind friendly compared with previous census
  figures.
- Add the rolling-rate atom-count heatmap as a candidate single-column
  companion figure. It uses a 10-year trailing detection-rate window and is
  intended to show time-local changes in discovery rate by atom-count class
  without relying on many overplotted traces. The fullerene class is labeled
  `Fuller` in the figure to avoid implying that the row is only C60.
- The detection-rate-by-atoms bubble figure has been changed from the 2021
  single-census view to a 2026-over-2021 comparison view. Explain that marker
  area now uses one absolute detections-to-area scale shared by both census
  views. This is a deliberate change from the 2021 legacy figure, where marker
  sizes were normalized within that single dataset and therefore were not
  suitable for direct cross-census comparison.
- The facility-share figure may be changed from the legacy pie grid to a
  horizontal bar chart. The bar chart is easier to compare quantitatively and
  labels each facility's denominator explicitly as `percent (n/N)`. If used,
  explain that each denominator is the number of first detections that occurred
  during that facility's operational lifetime in the selected census view.
- Swap the facility-share bar chart into the 2026 manuscript in place of the
  legacy pie-grid figure. The pie grid remains useful as an audit/reproduction
  artifact, but the bar chart should be the manuscript-facing version.
- The modern scopes-by-year figure emphasizes Yebes 40-m in MIT Red and uses
  faded dotted tails for NRAO 36-ft, Nobeyama 45-m, and NRAO 140-ft after their
  last credited first-detection contribution. If used, the caption should state
  that the dotted tails mean no later credited first detections in the selected
  census view, not necessarily physical telescope inactivity.
- A facility-era stacked bar chart was tested as a possible companion to the
  cumulative scopes-by-year trace and rejected for the 2026 manuscript.
- The periodic-table heatmap will retain the legacy yellow-to-red palette.
  `cividis` and `viridis` variants were previewed but rejected. The caption
  should specify that counts are secure, non-isotopologue ISM/CSM species and
  that lanthanides and actinides are omitted as in previous versions.
- Revisit the mass-by-wavelength figure caption after verification. The 2021
  caption said fullerenes were excluded, but the published IR sample count
  appears to include the three fullerene IR detections.
- Replace the 2021 mass-by-wavelength KDE in the 2026 manuscript with the
  horizontal box/whisker plus jittered-detection view. Keep the KDE generator in
  the package for historical comparison, but use the boxplot form because it
  better exposes the high-mass cm tail without implying smooth distributions
  for sparse IR, sub-mm, and UV/Vis samples. The production view should exclude
  fullerenes, matching the intended 2021 caption rather than the legacy plotting
  artifact. Caption should state that boxes show the interquartile range,
  whiskers span the 10th-90th percentiles, points are individual secure
  non-isotopologue non-fullerene ISM/CSM first detections, and the vertical
  marker highlights 80 amu.
- Replace the 2021 six-panel `mols_waves_by_atoms` KDE/histogram figure with
  the bubble heatmap form for the 2026 manuscript. Keep the six-panel generator
  for historical reproduction. The production bubble map uses the same
  yellow-to-red count palette as the periodic-table heatmap and should be
  described as counts of secure non-isotopologue non-fullerene ISM/CSM first
  detections by wavelength category and atom count. A row-percentage variant
  was considered but rejected in favor of absolute counts because the manuscript
  point is the observed population structure, not only the within-wavelength
  fractional distribution.
- Replace the 2021 DU histogram with the exact-value DU bar chart in the 2026
  manuscript. The histogram generator remains available for historical
  reproduction, but the bar chart better represents DU as a discrete/half-
  discrete formula-derived value rather than a continuous variable. The modern
  plot excludes negative-DU values from small protonated hydride/non-carbon
  species (`H2Cl+`, `H3+`, `H3O+`, and `NH3D+`), because those are artifacts of
  applying the neutral/organic DU expression outside its chemically meaningful
  domain. Caption should state that fullerenes are excluded and that only
  DU-compatible formulas containing H/D/N/C/O/S/F/Cl are included.
- Replace the 2021 DU-by-source-type KDE with the horizontal box/strip view in
  the 2026 manuscript. Keep the KDE generator for historical reproduction, but
  use the box/strip form because DU is a discrete/half-discrete formula-derived
  value and the manuscript point is better served by exact values plus
  medians/spreads than by smoothed overlapping KDEs. Caption should use the
  same summary language as the mass-by-wavelength boxplot: boxes show the
  interquartile range, whiskers span the 10th-90th percentiles, black lines
  mark medians, points are individual detections, fullerenes are excluded, and
  negative-DU domain artifacts are omitted from the production view. The 2026
  figure should replace the legacy KDE for the same manuscript point and uses
  Dark Cloud as the first row because it is now the largest high-DU source
  population.
- Replace the 2021 relative-DU-by-source-type KDE panel with the matching
  horizontal box/strip view in the 2026 manuscript. This companion plot should
  keep the same source-row order as the absolute DU-by-source figure and should
  be described as relative DU = DU/maxDU. Caption should note that fullerenes
  are excluded and negative-DU domain artifacts are omitted from the production
  view.
- Replace the 2021 mass-by-source-type KDE with the horizontal box/strip view
  in the 2026 manuscript. Keep the KDE generator for historical reproduction,
  but use the box/strip form because it shows exact molecular masses and
  summary spread without implying smooth source-specific mass distributions.
  Caption should state that each molecule is credited at most once per
  generalized first-detection source category, fullerenes are excluded, boxes
  show the interquartile range, whiskers span the 10th-90th percentiles, black
  lines mark medians, and points are individual secure non-isotopologue
  non-fullerene ISM/CSM first detections.
- Replace the 2021 wavelength-by-source-type pie grid with the normalized
  horizontal stacked-bar chart in the 2026 manuscript. The stacked bar preserves
  the same wavelength-credit data but makes cross-source comparison clearer.
  Caption should state that each molecule is credited at most once per
  generalized source category, every wavelength listed for that first detection
  is credited, bars are normalized within source category, and `n=` gives the
  number of wavelength credits, not necessarily the number of unique molecules.
- Implement the kappa histogram in the 2026 manuscript with the prolate to
  oblate guide arrow above the plot. Caption or surrounding prose should state
  that linear rotors are assigned the limiting value `kappa = -1`, nonlinear
  molecules require all three rotational constants, and molecules without a
  usable kappa are omitted from the plotted distribution.
- The molecule-type, source-type, and individual-source ring charts now sort
  rings dynamically by the selected census view's category counts while keeping
  stable category colors. For the 2026/current source-type chart, Dark Cloud
  has overtaken SFR and should appear as the outer ring. If comparing directly
  to 2021, note that ring order now reflects the current category ranking
  rather than a fixed historical order. The individual-source chart includes
  G+0.693-0.027 as `G+0.693` with a black inner ring and label; its labels are
  placed by rendered ring order so they stay aligned with the dynamically
  sorted rings. The individual-source chart counts first-detection source
  contributions, not unique molecules per source, so its percentages can sum
  to more than 100%.
- Replace the 2021 molecule-type-by-source pie grid with the enrichment-factor
  matrix in the 2026 manuscript. The matrix makes the same point more directly:
  source categories differ in which molecule classes are fractionally enriched
  or depleted relative to the overall secure ISM/CSM inventory. Caption should
  state that each cell shows an enrichment factor relative to the overall
  molecule-type fraction plus the raw count, and that molecule-type classes can
  overlap.

## Census Boundary Notes

- `history.accepted.census` is the authoritative field for whether a detection
  was accepted into a given census table. This distinction matters for species
  that were discussed as tentative or disputed in 2021 but are accepted in the
  2026 census.
- Until the 2026 census cutoff is frozen, the 2026 census view and the live
  current-database view are expected to match. After cutoff, post-2026 additions
  should remain available in current views but not in the frozen 2026 view.

## Associated Output Products

- Add a PPD detections slide generator in addition to the planned ISM molecule
  slide generator.
- Appendix/conceptual figures from the 2021 arXiv bundle were not migrated as
  database-generated products. Decide during manuscript drafting whether to
  reuse those static figures as-is or drop the appendix material entirely.
