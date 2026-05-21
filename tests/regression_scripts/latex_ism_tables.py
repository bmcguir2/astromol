import re
from pathlib import Path
from tempfile import TemporaryDirectory

from baseline import load_production_baseline

from astromol.census import CensusView
from astromol.database import Database
from astromol.latex import (
    BALANCED_ISM_MAX_COLUMNS,
    BALANCED_ISM_MAX_ROWS,
    balanced_ism_table_columns,
    ism_table_columns,
    ism_table_fragments,
    ism_table_molecules,
    write_ism_tables,
)


db = Database()
counts_2026 = load_production_baseline()["regression_counts"][
    "latex_ism_tables_2026"
]
view_2021 = CensusView.for_census(db, "2021")
view_2026 = CensusView.for_census(db, "2026")

two_seven, eight_more = ism_table_columns(view_2021)

assert [len(column) for column in two_seven] == [
    21,
    20,
    23,
    22,
    16,
    15,
    16,
    15,
    23,
    15,
]
assert [len(column) for column in eight_more] == [15, 14, 6, 6, 5, 2, 3, 3]
assert max(len(column) for column in two_seven) == 23
assert max(len(column) for column in eight_more) == 15

legacy_fragments = ism_table_fragments(view_2021, layout="legacy")
assert set(legacy_fragments) == {"ism_table_2-7.tex", "ism_table_8+.tex"}
assert legacy_fragments["ism_table_2-7.tex"].count(r"\molref{mol:") == 186
assert legacy_fragments["ism_table_8+.tex"].count(r"\molref{mol:") == 54

expected_first_two_seven = (
    r"\molref{mol:CH}{CH}"
    "\t&\t"
    r"\molref{mol:NH}{NH}"
    "\t&\t"
    r"\molref{mol:H2O}{H2O}"
    "\t&\t"
    r"\molref{mol:MgCN}{MgCN}"
    "\t&\t"
    r"\molref{mol:NH3}{NH3}"
    "\t&\t"
    r"\molref{mol:SiC3}{SiC3}"
    "\t&\t"
    r"\molref{mol:HC3N}{HC3N}"
    "\t&\t"
    r"\molref{mol:C4H-}{C4H-}"
    "\t&\t"
    r"\molref{mol:CH3OH}{CH3OH}"
    "\t&\t"
    r"\molref{mol:CH3CCH}{CH3CCH}\\"
)
expected_first_eight_more = (
    r"\molref{mol:HCOOCH3}{HCOOCH3}"
    "\t&\t"
    r"\molref{mol:CH3OCH3}{CH3OCH3}"
    "\t&\t"
    r"\molref{mol:CH3COCH3}{CH3COCH3}"
    "\t&\t"
    r"\molref{mol:HC9N}{HC9N}"
    "\t&\t"
    r"\molref{mol:C6H6}{C6H6}"
    "\t&\t"
    r"\molref{mol:C6H5CN}{C6H5CN}"
    "\t&\t"
    r"\molref{mol:1-C10H7CN}{1-C10H7CN}"
    "\t&\t"
    r"\molref{mol:C60}{C60}\\"
)

assert legacy_fragments["ism_table_2-7.tex"].splitlines()[8] == expected_first_two_seven
assert legacy_fragments["ism_table_8+.tex"].splitlines()[8] == expected_first_eight_more
assert r"\molref{mol:c-C3HCCH}{c-C3HCCH}" in legacy_fragments["ism_table_2-7.tex"]

linked_labels = set()
for content in legacy_fragments.values():
    linked_labels.update(re.findall(r"\\molref\{(mol:[^}]+)\}", content))

expected_labels = {molecule.label for molecule in view_2021.ism_molecules()}
assert linked_labels == expected_labels
assert len(linked_labels) == 240

with TemporaryDirectory() as tmp:
    output_dir = Path(tmp)
    assert write_ism_tables(view_2021, output_dir, layout="legacy") == legacy_fragments
    for filename, content in legacy_fragments.items():
        assert (output_dir / filename).read_text() == content

balanced_groups = balanced_ism_table_columns(view_2026)
assert len(balanced_groups) == 3
assert [len(group) for group in balanced_groups] == counts_2026[
    "balanced_group_column_counts"
]
assert all(len(group) <= BALANCED_ISM_MAX_COLUMNS for group in balanced_groups)
assert all(
    len(left) >= len(right)
    for left, right in zip(balanced_groups, balanced_groups[1:])
)
assert all(
    len(column.molecules) <= BALANCED_ISM_MAX_ROWS
    for group in balanced_groups
    for column in group
)

balanced_fragments = ism_table_fragments(view_2026)
assert set(balanced_fragments) == {
    "ism_table_1.tex",
    "ism_table_2.tex",
    "ism_table_3.tex",
}
assert (
    r"\multicolumn{2}{c}{\hyperref[2atoms]{2 Atoms}}"
    in balanced_fragments["ism_table_1.tex"]
)
assert (
    r"\multicolumn{3}{c}{\hyperref[3atoms]{3 Atoms}}"
    in balanced_fragments["ism_table_1.tex"]
)
assert (
    r"\multicolumn{2}{c}{\hyperref[5atoms]{5 Atoms}}"
    in balanced_fragments["ism_table_2.tex"]
)
assert r"\hyperref[13plusatoms]{13+ Atoms}" in balanced_fragments["ism_table_3.tex"]
assert r"\hyperref[14plusatoms]{14+ Atoms}" not in balanced_fragments["ism_table_3.tex"]

balanced_labels = []
for content in balanced_fragments.values():
    balanced_labels.extend(re.findall(r"\\molref\{(mol:[^}]+)\}", content))

expected_2026_molecules = ism_table_molecules(view_2026)
assert all(molecule.isotopologue_of is None for molecule in expected_2026_molecules)
expected_2026_labels = {molecule.label for molecule in expected_2026_molecules}
assert set(balanced_labels) == expected_2026_labels
assert len(balanced_labels) == len(expected_2026_labels) == counts_2026[
    "linked_labels"
]
assert len(balanced_labels) == len(set(balanced_labels))

with TemporaryDirectory() as tmp:
    output_dir = Path(tmp)
    assert write_ism_tables(view_2026, output_dir) == balanced_fragments
    for filename, content in balanced_fragments.items():
        assert (output_dir / filename).read_text() == content

print("LaTeX ISM table verification passed")
