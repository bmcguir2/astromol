# Data Model

Production data live in `astromol/data/`:

- `molecules.json`
- `detections.json`
- `sources.json`
- `telescopes.json`
- `references.bib`

`astromol.database.Database` loads these files, resolves cross-references, and
exposes the records as dataclass instances from `astromol.models`.

## Stable Identifiers

Records use stable identifier namespaces:

- molecules: `mol:<label>`
- detections: `det:<molecule>:<context>:<year>`
- sources: source `nick`
- telescopes: telescope `nick`
- references: Zotero/BibTeX citekeys

Detection IDs are stable links used for disputed, tentative, confirming, and
confirmed-by relationships.

## History Metadata

Molecules, detections, sources, and telescopes carry record-history metadata.
The history distinguishes:

- when a record was introduced into tracking
- when it was accepted into a census inventory
- when it was last modified

This allows historical 2018 and 2021 views, the developing 2026 census view, and
the live current database view to be computed from the same production records.

## Census Views

`astromol.census.CensusView` centralizes scientific selection rules for output
generation:

- `CensusView.for_census(db, "2021")` freezes records at a census boundary.
- `CensusView.for_census(db, "2026")` represents the developing 2026 census.
- `CensusView.current(db)` represents the live database.

Secure detections are selected by accepted history. Tentative and disputed
detections are selected only when explicitly requested.

Isotopologues are excluded from context views by default. Pass
`include_isotopologues=True` when isotope-expanded inventories are needed.
