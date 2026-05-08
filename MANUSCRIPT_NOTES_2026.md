# 2026 Manuscript Notes

This file collects writing-time reminders that arose during the database and
output-generation migration. Keep implementation details in `SPEC.md` and
`GENERATION_MIGRATION.md`; keep manuscript-facing explanations here.

## Terminology And Classification

- Explain the source-type terminology change from `LOS Cloud` to
  `Diffuse Cloud`. The new term is meant to describe diffuse-cloud
  line-of-sight material more directly, while still allowing the source table to
  consolidate these sight lines into a single `Diffuse Cloud` row.

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
