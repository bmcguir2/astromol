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

### Molecule

Intrinsic chemical species metadata.

Required fields:
- `name`
- `formula`

Display and notes:
- `table_formula` (defaults to `formula`)
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

Spectroscopy:
- `rotcon` (`RotationalConstants` or null)
- `dipole` (`DipoleMoment` or null)

References and relationships:
- `refs`: role-keyed dict using roles from `MOLECULE_REF_ROLES`
- `isotopologue_of`: parent molecule label or null

LaTeX fields:
- `latex_header`
- `latex_notes`

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
- `molecule`: molecule label, resolved to a `Molecule`
- `sources`: source nicks, resolved to `Source` objects
- `telescopes`: telescope nicks, resolved to `Telescope` objects
- `wavelengths`: values from `WAVELENGTHS`
- `year`
- `type`: value from `DETECTION_TYPES`

Optional fields:
- `note`
- `first`
- `month`
- `day`
- `refs`: role-keyed dict using roles from `DETECTION_REF_ROLES`
- `latex_text`

Computed:
- `sortdate`

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
