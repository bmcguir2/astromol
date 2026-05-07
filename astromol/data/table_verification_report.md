# 2021 Table Verification Report

This temporary audit report records cross-verification of the refactored
`astromol` JSON data against LaTeX tables from the 2021 census manuscript.

The purpose is to verify that the refactor preserved the curated 2021 census
content while also documenting places where the modern schema intentionally
differs from the historical table representation.

## Verification Rules

- Source tables are the `.tex` files staged in `astromol/data/`.
- Database checks use current production JSON files unless otherwise noted.
- "2021 accepted census population" means molecules with
  `history.accepted.census <= 2021`.
- For known ISM/CSM molecule tables, database-side detections are filtered to
  confirmed `type: "ISM/CSM"` detections.
- Tentative and disputed detections are handled according to the modern
  `detections.json` status model. Their molecule records may still exist in the
  database, but they are not treated as accepted census entries unless
  `history.accepted` records an accepted census.
- Context-specific tables, such as ice or protoplanetary-disk tables, are
  checked against detections in that context. Secure context-specific table
  membership is selected with detection `history.accepted.census`; detection
  year alone is not sufficient because census cutoffs can fall partway through
  a calendar year. When a historical table includes a tentative/disputed item,
  the audit reports both exact historical reproduction and the modern
  secure-only interpretation.
- A complete 2018-vs-2021 accepted-history reconstruction is deferred until
  the 2018 census tables are added to this audit. Detection-level
  `history.accepted.census` is reliable for contexts already audited here, but
  should not be treated as globally authoritative for the 2018/2021 boundary
  until the corresponding 2018 tables have been cross-checked.
- Differences are classified before any corrective action:
  - `conversion_error`: current data failed to preserve intended legacy content.
  - `intentional_schema_change`: modern identifiers or fields intentionally
    differ from the 2021 table presentation.
  - `status_formalism_change`: modern detection status or accepted-census
    modeling changes table membership relative to the historical manuscript
    table.
  - `legacy_table_issue`: the 2021 table itself appears incomplete, inconsistent,
    or superseded.
  - `needs_human_decision`: discrepancy needs curator judgment.

## ISM Tables

Input files:
- `astromol/data/ism_table_2-7.tex`
- `astromol/data/ism_table_8+.tex`

Check:
- Parsed every `\hyperref[...]{\ce{...}}` molecule entry from both tables.
- Mapped legacy/custom table labels to current molecule labels where already
  decided during refactor, e.g. `CNN1 -> mol:1-C10H7CN`,
  `CNN2 -> mol:2-C10H7CN`, and `SiC2 -> mol:c-SiC2`.
- Compared against confirmed current `ISM/CSM` detections restricted to
  molecules accepted in or before the 2021 census according to
  `history.accepted.census`.
- Checked table category assignment against current atom counts, PAH flags, and
  fullerene flags.

Results:
- Parsed expected table entries: 240.
- Expected entries mapped to current molecule labels: 240/240.
- Current confirmed `ISM/CSM` entries with `history.accepted.census <= 2021`:
  240.
- Missing expected entries in current confirmed 2021-scope ISM/CSM set: 0.
- Extra current confirmed 2021-scope ISM/CSM entries absent from the historical
  tables: 0.
- Atom-count/category mismatches: 0.
- PAH/fullerene classification mismatches: 0.

Reproduced table category counts:

| Category | Count |
| --- | ---: |
| 2 atoms | 41 |
| 3 atoms | 45 |
| 4 atoms | 31 |
| 5 atoms | 31 |
| 6 atoms | 23 |
| 7 atoms | 15 |
| 8 atoms | 15 |
| 9 atoms | 14 |
| 10 atoms | 6 |
| 11 atoms | 6 |
| 12 atoms | 5 |
| 13 atoms | 2 |
| PAHs | 3 |
| Fullerenes | 3 |

Discrepancies and notes:

| Item | Historical Table | Current Data | Classification | Resolution |
| --- | --- | --- | --- | --- |
| Silacyclopropynylidene | `SiC2` | `mol:c-SiC2`, `table_formula: c-SiC2` | `intentional_schema_change` | Accepted. The cyclic prefix was restored during refactor because it is chemically meaningful. |

Previously flagged status/formalism warnings now resolved by
`history.accepted`:

| Item | `history.introduced` | `history.accepted` | Resolution |
| --- | --- | --- | --- |
| Ethyl methyl ether (`mol:C2H5OCH3`) | `census: 2021`, `context: tentative` | `census: 2026`, `context: confirmed_ism_csm` | Excluded from historical 2021 accepted table while retained in the database with its later confirmed detection. |
| N-methyl formamide (`mol:CH3NHCHO`) | `census: 2021`, `context: tentative` | `census: 2026`, `context: confirmed_ism_csm` | Excluded from historical 2021 accepted table while retained in the database with its later confirmed detection. |
| 1-butyne (`mol:CH3CH2CCH`) | `census: 2021`, `context: tentative` | `census: 2026`, `context: confirmed_ism_csm` | Excluded from historical 2021 accepted table while retained in the database with its later confirmed detection. |

Conclusion:

The two historical known-ISM tables are fully reproduced by the current data
when interpreted as historical 2021 accepted known-molecule tables. The only
remaining difference is the intentional `SiC2` to `c-SiC2` schema correction,
not a conversion error.

## Ice Table

Input file:
- `astromol/data/ice_table.tex`

Check:
- Parsed every `\ce{...}` species entry from the table body.
- Mapped table formulas to current molecule labels.
- Compared against current `type: "ice"` detections dated in or before 2021.
- Checked the historical table reproduction both with all statuses and with
  modern secure-only filtering.
- Checked that representative observation references are present in
  `references.bib` under the normalized citekeys used by `detections.json`.

Results:
- Parsed expected table entries: 10.
- Expected entries mapped to current molecule labels: 10/10.
- Current ice detections dated in or before 2021, all statuses: 10.
- Missing expected entries in the all-status historical reproduction: 0.
- Extra current entries in the all-status historical reproduction: 0.
- Current secure ice detections dated in or before 2021: 9.
- Missing expected entries in the modern secure-only interpretation: 1
  (`mol:OCN-`).
- Extra current entries in the modern secure-only interpretation: 0.
- Current post-2021 secure ice detections excluded from the historical table:
  5 (`mol:CH3CHO`, `mol:HCOOCH3`, `mol:CH3OCH3`, `mol:CH3CH2OH`,
  `mol:CH3COCH3`).

Reference mapping:

| Historical key | Current key |
| --- | --- |
| `1979ApJ...232L..53S` | `Soifer:1979:L53` |
| `1973ApJ...179..483G` | `Gillett:1973:483` |
| `1995ApJ...449..674P` | `Palumbo:1995:674` |
| `1997ApJ...479..839P` | `Palumbo:1997:839` |
| `1989A&A...223L...5D` | `DHendecourt:1989:L5` |
| `2005A&A...441..249V` | `vanBroekhuizen:2005:249` |
| `1998ApJ...501L.105L` | `Lacy:1998:L105` |
| `2001A&A...376..254K` | `Keane:2001:254` |
| `1999A&A...343..966S` | `Schutte:1999:966` |
| `1991ApJ...376..556L` | `Lacy:1991:556` |
| `1991A&A...243..473G` | `Grim:1991:473` |

Discrepancies and notes:

| Item | Historical Table | Current Data | Classification | Resolution |
| --- | --- | --- | --- | --- |
| Cyanate anion | `OCN-` included in the ice table | `mol:OCN-`, `det:OCN-:ice:2005`, `status: tentative`, `history.accepted: null` | `status_formalism_change` | Curator accepted this difference as a 2021 simplicity choice. Historical reproduction includes it because the 2021 table printed it; 2026 generation will handle the tentative status in prose rather than changing the database. |

Conclusion:

The historical 2021 ice table is fully reproduced by the current data when the
audit includes the same tentative item printed in that table. Under the modern
secure-only detection model, `OCN-` is intentionally excluded; this is a
status/formalism distinction, not a conversion error.

## PPD Table

Input file:
- `astromol/data/ppd_table.tex`

Check:
- Parsed every `\ce{...}` species entry from the table body.
- Mapped table formulas, including mhchem isotope notation, to current molecule
  labels via `table_formula`.
- Compared against current secure `type: "ppd"` detections accepted in or
  before the 2021 census according to detection `history.accepted.census`.

Results:
- Parsed expected table entries: 40.
- Expected entries mapped to current molecule labels: 40/40.
- Current secure PPD detections with `history.accepted.census <= 2021`: 40.
- Missing expected entries: 0.
- Extra current entries: 0.
- Current post-2021 secure PPD detections excluded from the historical table:
  12 (`mol:13CO2`, `mol:33SO`, `mol:34SO`, `mol:C6H6`, `mol:CH3OCH3`,
  `mol:H13CCH`, `mol:H213CO`, `mol:HC18O+`, `mol:HC4H`, `mol:NO`,
  `mol:SiS`, `mol:c-C2H4O`).

Modeling issue resolved by this table:

The PPD table is context-specific. Molecule-level `history.accepted` is not
sufficient to reproduce context-specific census tables because a molecule can be
accepted in the ISM while its PPD detection is not yet part of a printed census
table. Conversely, several isotopologue molecule records printed in the 2021
PPD table were introduced as standalone records during the refactor and still
carry molecule `history.accepted.census: 2026`, even though their PPD detection
membership is historically 2021. This was resolved by extending
`history.accepted` to detections and using detection-level accepted census
membership for context-specific table generation.

Discrepancies and notes:

| Item | Historical Table | Current Data | Classification | Resolution |
| --- | --- | --- | --- | --- |
| PPD isotopologues printed in 2021 | `C^{15}N`, `DCO+`, `H^{13}CO+`, `^{13}CO`, `DCN`, `C^{18}O`, `H^{13}CN`, `C^{17}O`, `H^{15}CN`, `HD`, `DNC`, `C^{34}S`, `^{13}CS`, `N2D+`, `C2D` | Molecule records may have `history.accepted.census: 2026` because they were promoted from nested/non-standalone legacy records during refactor; their PPD detections now have `history.accepted.census: 2021`. | `status_formalism_change` | Resolved with detection-level accepted census membership. The PPD table is generated from detection history, not molecule history. |
| Secure PPD detections absent from 2021 table | `13C17O`, `13C18O`, `C2S`, `CH2CN`, `SO2` absent | Current detections are secure PPD detections with years `2017`, `2019`, or `2021`; all now have detection `history.accepted.census: 2026`. | `status_formalism_change` | Curator confirmed these were added to the database after the 2021 census cutoff. Detection year alone is insufficient; detection accepted census records the correct table membership. |

Conclusion:

The historical 2021 PPD table is fully reproduced by selecting secure PPD
detections with `history.accepted.census <= 2021`. Detection-level accepted
census membership resolves both the isotopologue standalone-record formalism
and the 2021 cutoff cases.

## Extragalactic Table

Input file:
- `astromol/data/exgal_table.tex`

Check:
- Parsed every `\ce{...}` species entry from the table body, including dagger
  markers for tentative detections.
- Mapped table formulas to current molecule labels via `table_formula`.
- Compared against current secure `type: "exgal"` detections accepted in or
  before the 2021 census according to detection `history.accepted.census`.
- Compared the full historical table, including tentative dagger-marked rows,
  against secure accepted entries plus tentative entries introduced by the 2021
  census.

Results:
- Parsed expected table entries: 75.
- Expected entries mapped to current molecule labels: 75/75.
- Historical table entries marked tentative with dagger symbols: 2
  (`mol:c-C3H`, `mol:HC5N`).
- Current secure exgal detections with `history.accepted.census <= 2021`: 73.
- Missing expected secure entries: 0.
- Extra current secure entries: 0.
- Missing expected entries in the full historical all-status reproduction: 0.
- Extra current entries in the full historical all-status reproduction: 0.

Discrepancies and notes:

| Item | Historical Table | Current Data | Classification | Resolution |
| --- | --- | --- | --- | --- |
| Post-cutoff 2021 exgal detections | `CH3CH2OH`, `H2NC`, `HC2CHO`, and `HOCN` absent | Current secure exgal detections have `year: 2021`; curator confirmed they entered after the 2021 census cutoff. | `status_formalism_change` | Corrected. Detection `history.accepted.census` is now `2026` for `det:CH3CH2OH:exgal:2021`, `det:H2NC:exgal:2021`, `det:HC2CHO:exgal:2021`, and `det:HOCN:exgal:2021`. |
| Tentative exgal entries | `c-C3H` and `HC5N` marked with dagger symbols | The legacy conversion preserved `note: Legacy detection flag: Tentative`, but the first detection-history normalization treated missing `status` as secure and marked both as accepted in 2021. | `conversion_error` | Corrected. `det:c-C3H:exgal:2006` and `det:HC5N:exgal:2015` now have `status: tentative`, `first: false`, `history.introduced.census: 2018`, `history.introduced.context: tentative`, and `history.accepted: null`. |

Cause of the 2021-vs-2018 accepted-history error:

Detection-level history was generated after the molecule-level history pass.
The temporary bulk rule used for contexts other than ISM/CSM, PPD, and ice
treated secure detections with `year <= 2021` as accepted in the 2021 census.
Because these two exgal records did not have structured `status: tentative`
yet, they were interpreted as secure and assigned `history.accepted.census:
2021`. The legacy tentative information existed only as free text in `note`,
so it was not used by the bulk history normalization.

Conclusion:

The extragalactic table is fully reproduced. Secure accepted 2021-scope exgal
detections reproduce the secure historical table entries, and the full
historical all-status reproduction also matches once the two dagger-marked
tentative detections are represented with structured status.

## Exoplanet Table

Input file:
- `astromol/data/exo_table.tex`

Check:
- Parsed every `\ce{...}` species entry from the table body.
- Mapped table formulas to current molecule labels via `table_formula`.
- Compared against current secure `type: "exo"` detections accepted in or
  before the 2021 census according to detection `history.accepted.census`.
- The table caption explicitly excludes tentative and disputed detections, so
  the database-side comparison is secure-only.

Results:
- Parsed expected table entries: 9.
- Expected entries mapped to current molecule labels: 9/9.
- Current secure exoplanet detections with `history.accepted.census <= 2021`:
  9.
- Missing expected entries: 0.
- Extra current secure entries: 0.

Discrepancies and notes:

| Item | Historical Table | Current Data | Classification | Resolution |
| --- | --- | --- | --- | --- |
| `13CO` exoplanet detection | Absent from the 2021 exoplanet table | `det:13CO:exo:2021` has `year: 2021` and `status: secure`; curator confirmed it entered after the 2021 census cutoff. | `status_formalism_change` | Corrected. Detection `history.accepted.census` is now `2026`. |

Conclusion:

The exoplanet table is fully reproduced by selecting secure exoplanet
detections accepted in or before the 2021 census.

## Source Count Table

Input file:
- `astromol/data/source_table.tex`

Check:
- Parsed every source/count pair from the two-column table body.
- Counted source contributions from current secure `type: "ISM/CSM"`
  detections accepted in or before the 2021 census according to detection
  `history.accepted.census`.
- Used current source `latex_name` values for display-name matching.
- Applied the historical table grouping rule that all diffuse-cloud/line-of-
  sight detections are consolidated as `LOS Cloud`. This is a historical table
  label; the modern source type is `Diffuse Cloud`.

Results:
- Parsed expected source rows: 38.
- Current grouped source rows: 38.
- Missing expected rows: 0.
- Extra current rows: 0.
- Count mismatches: 0.

Conclusion:

The source contribution table is fully reproduced by the current data when the
historical `LOS Cloud` grouping convention is applied to modern
`Diffuse Cloud` source records.

## Facilities Count Table

Input file:
- `astromol/data/facilities_table.tex`

Check:
- Parsed every facility/count pair from the two-column table body.
- Counted facility contributions from current secure `type: "ISM/CSM"`
  detections accepted in or before the 2021 census according to detection
  `history.accepted.census`.
- Used current telescope `latex_name` values for display-name matching.

Results:
- Parsed expected facility rows: 46.
- Current facility rows: 46.
- Missing expected rows: 0.
- Extra current rows: 0.
- Count mismatches: 0.

Conclusion:

The facility contribution table is fully reproduced by the current data.
