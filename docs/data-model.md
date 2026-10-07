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

Each record has a stable identifier so other records can link to it:

- molecules: `mol:<label>`
- detections: `det:<molecule>:<context>:<year>`
- sources: source `nick`
- telescopes: telescope `nick`
- references: Zotero/BibTeX citekeys

Detection relationships use these IDs to link claims, confirmations, disputes,
and superseding records.

## History Metadata

Molecules, detections, sources, and telescopes carry record-history metadata.
The history distinguishes:

- when a record was introduced into tracking
- when it was accepted into a census inventory
- when it was last modified

This allows published 2018 and 2021 membership views and the live working
inventory for the next census to be computed from the same production records.

## Census Views

`astromol.census.CensusView` applies the same membership and filtering rules to
each output:

- `CensusView.for_census(db, "2021")` selects membership at a published census boundary.
- `CensusView.current(db)` selects the live working inventory.

Secure detections are selected by accepted history. Tentative and disputed
detections are selected only when explicitly requested.

Context views exclude isotopologues by default. Pass
`include_isotopologues=True` to include them.
