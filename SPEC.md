# astromol Architecture Specification

## Overview

astromol is a Python library for cataloging astronomical molecule detections.
Data is stored in JSON files, with bibliographic metadata loaded from a Zotero
BibTeX export. The code loads data, resolves cross-references, and provides
query helpers through `Database`.

## File Structure

```text
astromol/
├── data/
│   ├── molecules.json
│   ├── detections.json
│   ├── sources.json
│   ├── telescopes.json
│   └── references.bib
├── models.py
└── database.py
```

`test_load.py` is a local smoke script. The data files may be incomplete during
curation, but loaded cross-references must resolve for `Database()` to succeed.

## Data Provenance

The current `molecules.json`, `detections.json`, and `sources.json` baseline was
promoted from the legacy conversion workflow in commit `67f37d9`. The final
audited preview state is preserved in commit `0187018`; use those commits for
the removed conversion scripts, staging reports, and legacy source data.

Molecule history records were imported from the legacy `census_version` and
`change_log` fields through a staged preview workflow. Legacy versions are
preserved as `history.events[*].legacy_version`, while
`history.introduced.census` records the first print census appearance. Molecules
present in the legacy database but absent from `2021_census_arxiv.tex` are
therefore marked as introduced in the `2026` census.

OpenAI Codex is used as a coding assistant for software implementation,
refactoring, validation scripts, documentation drafts, and workflow support.
Scientific database content is not AI-generated. Molecule records, detection
records, source/telescope metadata, reference mappings, classifications, and
curation decisions are derived from the scientific literature, Zotero-managed
bibliography exports, legacy census materials, and human curator review. AI
assistance may be used to stage or transform records for review, but no staged
scientific data are accepted into production without human verification.

## Curation Workflow

Production data lives in the JSON files under `astromol/data/`. New records
should normally be staged first using copy/paste YAML templates from
`curation/templates/`, with working staging files kept under `curation/staging/`.
Run `python scripts/stage_records.py --staging <file.yaml>` to generate preview
JSON and a review report. Use `--apply` only after the preview is accepted.
The staging script automatically populates semantic history metadata for new
records, including `history.introduced.date`, `history.last_modified`, and an
initial `added` event dated with the staging run date.

When the schema for `Molecule`, `Detection`, `Source`, or `Telescope` changes,
the corresponding curation template, `scripts/stage_records.py`, and this
specification must be updated in the same change.

## Data Model

### Ref

Bibliographic reference loaded from `references.bib`.

Fields:
- `bibcode` (str, required): BibTeX citation key using
  `LastName:Year:FirstPage`, e.g. `Swings:1937:483`. Add lowercase suffixes
  only when needed to disambiguate duplicate keys, e.g. `Cabezas:2021:L9a`.
- `author` (str, required): compact author display string
- `journal` (str, required): journal, booktitle, publisher, or entry type
- `year` (int, required)
- `month` (int or null)
- `day` (int or null)
- `volume` (str or null)
- `page` (str or null)
- `title` (str or null)
- `doi` (str or null)
- `note` (str or null): currently populated from BibTeX `annotation`

Computed:
- `sortdate`

### Telescope

Observing facility or instrument.

Required fields:
- `name`
- `nick`
- `shortname`
- `type`
- `wavelength` (list of values from `WAVELENGTHS`)

Optional fields:
- `diameter`
- `latitude`
- `longitude`
- `built`
- `decommissioned`
- `note`
- `latex_name`
- `history`: optional `RecordHistory`

Computed:
- `active`

### Source

Astronomical object where molecules have been detected.

Required fields:
- `name`
- `nick`
- `type` (one of `SOURCE_TYPES`)

Optional fields:
- `ra`
- `dec`
- `simbad_url`
- `latex_name`
- `note`
- `history`: optional `RecordHistory`

Current `SOURCE_TYPES`:
- `SFR`
- `Dark Cloud`
- `Carbon Star`
- `Oxygen Star`
- `Planetary Nebula`
- `Protostellar`
- `Diffuse Cloud`
- `HII Region`
- `PDR`
- `Shock`
- `Supernova Remnant`
- `Sgr A`
- `External Galaxy`
- `YSO`
- `Exoplanet`
- `Other`

### RotationalConstants

Rotational constants in MHz.

Fields:
- `A`
- `B`
- `C`
- `refs` (list): BibTeX keys in JSON, resolved to `Ref` objects during loading
- `note`

Compatibility:
- `ref` returns the first item in `refs`, or `None`.

### DipoleMoment

Dipole moment components in Debye.

Fields:
- `a`
- `b`
- `c`
- `refs` (list): BibTeX keys in JSON, resolved to `Ref` objects during loading
- `note`

Computed:
- `total`

Compatibility:
- `ref` returns the first item in `refs`, or `None`.

### RecordHistory

Semantic update history exposed through database records. This is not intended
to replace git as the fine-grained audit trail. It records user-facing database
events that explain when and why a record was added, updated, or corrected.

History fields:
- `introduced`: dict with optional `date` and `census` keys
- `last_modified`: ISO date string for the last semantic database change
- `events`: list of `HistoryEvent`

`HistoryEvent` fields:
- `kind`: one of `added`, `updated`, or `corrected`
- `summary`: concise human-readable description of the event
- `date`: optional ISO date string
- `census`: optional census release string, such as `2021` or `2026`
- `fields`: optional list of affected field names
- `legacy_version`: optional legacy astromol version identifier preserved from
  the pre-refactor database

Records introduced in a given census can be selected with
`history.introduced.census`. Records changed for a given census can be selected
from `history.events[*].census`, independently of whether they were newly
introduced.

### Molecule

Intrinsic chemical species metadata.

Required fields:
- `name`
- `formula`: canonical internal chemical formula. This must be parseable by
  `molmass`; isotopologues use explicit bracket notation such as `H2[13C]O`.

Display and notes:
- `table_formula`: readable formula body for manuscript tables. This may use
  mhchem isotope placement such as `H2^{13}CO`, but it does not include the
  surrounding `\ce{}` wrapper; LaTeX generation adds that wrapper.
  Defaults to `formula`.
- `label` (required by `__post_init__`)
- `note`

Optional identifiers:
- `iupac_name`
- `selfies`
- `synonyms`
- `smiles`
- `canonical_smiles`
- `inchi`
- `inchikey`

Curated flags and structure metadata:
- `radical_override`
- `fullerene`
- `pah`
- `n_rings`
- `cyclic`
- `tags`: dictionary of manually curated analysis annotations. Values are lists
  of strings. Current category:
  - `functional_groups`: manually assigned functional-group or structural-motif
    tags used for analysis. These are curated annotations, not computed
    structural assertions.

Functional-group tag definitions:
- `CN`: Cyano/nitrile-family C-N motif. Includes nitriles, metal cyanides,
  cyanopolyynes, and related species manually classified as containing a
  cyano-like C-N unit.
- `NC`: Isocyano/isocyanide-family N-C motif, where the connectivity is the
  structural isomer of a cyano/nitrile-like group.
- `CH3`: Methyl group or methyl substituent.
- `=O`: Carbonyl/oxygen multiple-bond motif, broadly including aldehydes,
  ketones, acids, esters, amides, ureas, ketenes, isocyanates, and
  O-terminated carbon-chain analogues.
- `=S`: Thiocarbonyl/sulfur multiple-bond motif, broadly the sulfur analogue
  of `=O`, including thioaldehydes, thioketones, thioketenes,
  isothiocyanates, and S-terminated carbon-chain analogues.
- `OH`: Hydroxyl group, including alcohols, hydroxy-substituted species,
  carboxylic acids, and metal hydroxides.
- `SH`: Thiol/mercapto group, including sulfur analogues of hydroxyl-bearing
  species.
- `NH2`: Amino group or primary amide/amine-like `NH2` substituent.
- `NH`: Imine/imino/secondary amine-like `NH` motif, manually assigned where
  `NH` is structurally meaningful rather than just present in the formula.

Spectroscopy:
- `rotcon` (`RotationalConstants` or null)
- `dipole` (`DipoleMoment` or null)

References and relationships:
- `refs`: role-keyed dict using roles from `MOLECULE_REF_ROLES`
  - `lab` includes laboratory characterization references for the molecule.
    Legacy context-specific lab references, such as ice laboratory spectra, are
    normalized here rather than attached to detections.
- `isotopologue_of`: parent molecule label or null

LaTeX fields:
- `latex_section_override`: optional custom section heading. This should be
  null for normal generated section headings.
- `latex_body`: curated prose for the molecule's manuscript section.

History:
- `history`: optional `RecordHistory`

Computed from formula or spectroscopy:
- `atom_counts`
- `atoms`
- `isotope_counts`
- `mass`
- `average_mass`
- `nominal_mass`
- `natoms`
- `charge`
- `cation`
- `anion`
- `neutral`
- `nelectrons`
- `odd_electron`
- `radical`
- `du`
- `maxdu`
- `is_linear`
- `kappa`

When a formula uses an isotope label that is valid for the database but
unsupported by `molmass`, the model uses `ISOTOPE_MASS_OVERRIDES` for hardcoded
relative atomic masses. Current override:
- `26Al`: 25.98689188 u

Loader compatibility:
- Legacy stored/computed fields such as `mass`, `natoms`, `charge`, `du`,
  `maxdu`, `is_linear`, and `kappa` are ignored before constructing
  `Molecule`.
- Legacy `radical` is mapped to `radical_override`.
- Legacy flat spectroscopy keys `Acon/Bcon/Ccon` and `mua/mub/muc` are mapped
  into `rotcon` and `dipole`.

### Detection

Specific molecule detection event.

Required fields:
- `id`: stable detection identifier
- `molecule`: molecule label, resolved to a `Molecule`
- `sources`: source nicks, resolved to `Source` objects
- `telescopes`: telescope nicks, resolved to `Telescope` objects
- `wavelengths`: values from `WAVELENGTHS`
- `year`
- `type`: value from `DETECTION_TYPES`

Optional fields:
- `note`
- `status`: value from `DETECTION_STATUSES`; defaults to `secure`
- `status_note`: concise explanation for tentative or disputed detections
- `first`
- `month`
- `day`
- `refs`: role-keyed dict using roles from `DETECTION_REF_ROLES`
- `confirms`: detection IDs this record confirms
- `confirmed_by`: detection IDs that confirm this record
- `disputes`: detection IDs this record disputes
- `disputed_by`: detection IDs that dispute this record
- `supersedes`: detection IDs this record supersedes
- `superseded_by`: detection IDs that supersede this record
- `latex_text`
- `history`: optional `RecordHistory`

Computed:
- `sortdate`

Detection IDs are stable database identifiers. The default format is:

```text
det:<molecule-label-without-mol-prefix>:<context>:<year>[:qualifier]
```

Examples:
- `det:CH:ism-csm:1937`
- `det:CH3CH2CCH:ism-csm:2021`
- `det:CH3CH2CCH:ism-csm:2024`
- `det:OCN-:ice:2005`

Use an additional qualifier only when molecule, context, and year would not be
unique. Detection IDs should not include mutable status terms such as
`tentative`, because a detection claim may later be confirmed without becoming a
different database record.

Detection status is independent of detection type. `type` records the
astrophysical context, such as `ISM/CSM`, `ice`, or `exgal`; `status` records
whether the claim is `secure`, `tentative`, or `disputed`. Molecule inclusion
does not imply a secure astronomical detection. Tentative and disputed records
should normally keep `first: false` unless a deliberate curation decision says
otherwise.

`first: true` means the first accepted/secure detection in that context, not
necessarily the first chronological claim. Earlier tentative or disputed claims
should be represented as separate detection records and linked to later secure
detections with `confirmed_by`/`confirms` where appropriate.

Relationship fields store detection IDs, not BibTeX keys. Use `refs` for papers
and relationship fields for links between detection records. For example, an
earlier tentative butyne detection can use
`confirmed_by: [det:CH3CH2CCH:ism-csm:2024]`, while the confirming paper remains
under `refs.confirmation` or the secure detection's `refs.observation`.

Current `DETECTION_TYPES`:
- `ISM/CSM`
- `isotopologue`
- `ice`
- `exgal`
- `ppd`
- `exo`
- `comet`

Current `DETECTION_STATUSES`:
- `secure`
- `tentative`
- `disputed`

Current `DETECTION_REF_ROLES`:
- `observation`: paper reporting the detection claim
- `confirmation`: paper confirming an earlier tentative detection
- `dispute`: paper challenging a detection claim
- `correction`: paper correcting detection metadata or interpretation

Current detection relationship fields:
- `confirms`
- `confirmed_by`
- `disputes`
- `disputed_by`
- `supersedes`
- `superseded_by`

## Reference Loading

`Database._load_refs()` reads `astromol/data/references.bib` using
`bibtexparser` and stores `Ref` objects in `Database.refs`, keyed by the BibTeX
citation key.

JSON references should use those BibTeX keys. The loader resolves:
- molecule-level `refs`
- detection-level `refs`
- nested spectroscopy `ref` or `refs`

Nested spectroscopy supports both current `refs: [...]` and legacy `ref`.

## Database

`Database()` loads data in this order:
1. references
2. telescopes
3. sources
4. molecules
5. detections

Primary stores:
- `refs`: dict keyed by BibTeX citation key
- `telescopes`: dict keyed by telescope nick
- `sources`: dict keyed by source nick
- `molecules`: dict keyed by unique molecule label
- `molecules_by_formula`: dict keyed by chemical formula, with list values
- `molecules_by_name`: dict keyed by molecule name, with list values
- `detections`: list of `Detection`

Molecule labels are the database identity keys. Labels should use the `mol:`
namespace, e.g. `mol:CH` or `mol:l-C3H`. Chemical formulae are data fields and
are not required to be unique because structural isomers can share a formula.
Detection `molecule` values and `isotopologue_of` values refer to molecule
labels, not formulae.

Accessors:
- `get_ref`
- `get_telescope`
- `get_source`
- `get_molecule`
- `get_molecules_by_formula`
- `get_molecules_by_name`

`Database.__repr__()` prints a compact count summary of loaded molecules,
detections, sources, telescopes, and references.

## Validation Expectations

For a complete dataset:
- all JSON files parse
- required dataclass fields are present
- all source types, detection types, and wavelengths are allowed
- molecule, source, telescope, and reference cross-references resolve
- molecule labels, source nicks, telescope nicks, and BibTeX keys are unique
- all `refs` entries use valid role names
