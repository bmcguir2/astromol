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
