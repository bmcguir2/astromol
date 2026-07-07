import re
from pathlib import Path
from tempfile import TemporaryDirectory

from baseline import load_production_baseline

from astromol.census import CensusView
from astromol.database import Database
from astromol.latex import (
    ppd_table_atom_groups,
    ppd_table_columns,
    ppd_table_detections,
    ppd_table_fragments,
    write_ppd_table,
)


db = Database()
counts_2026 = load_production_baseline()["regression_counts"]["latex_ppd_table_2026"]
view_2021 = CensusView.for_census(db, "2021")
view_2026 = CensusView.for_census(db, "2026")

assert ppd_table_atom_groups(view_2021) == [(2, 3, 4, 5, 6)]

columns_2021 = ppd_table_columns(view_2021)
assert [[len(column) for column in group] for group in columns_2021] == [
    [14, 16, 4, 4, 2],
]

detections_2021 = ppd_table_detections(view_2021)
assert len(detections_2021) == 40
assert sum(detection.molecule.isotopologue_of is not None for detection in detections_2021) == 15
assert all(detection.status == "secure" for detection in detections_2021)
assert all(detection.refs.get("observation") for detection in detections_2021)

fragments_2021 = ppd_table_fragments(view_2021)
assert set(fragments_2021) == {"ppd_table.tex"}
content_2021 = fragments_2021["ppd_table.tex"]
assert content_2021.count(r"\molref{mol:") == 40
assert content_2021.count(r"\citet{") == 38
assert r"$^{\dagger}$" not in content_2021

expected_first_row = (
    r"\molref{mol:CN}{CN}"
    "\t&\t1, 2\t&\t"
    r"\molref{mol:H2O}{H2O}"
    "\t&\t3, 4, 5\t&\t"
    r"\molref{mol:NH3}{NH3}"
    "\t&\t6\t&\t"
    r"\molref{mol:HC3N}{HC3N}"
    "\t&\t7\t&\t"
    r"\molref{mol:CH3OH}{CH3OH}"
    "\t&\t8"
    r"\\"
)
assert content_2021.splitlines()[9] == expected_first_row

linked_labels_2021 = set(re.findall(r"\\molref\{(mol:[^}]+)\}", content_2021))
expected_labels_2021 = {detection.molecule.label for detection in detections_2021}
assert linked_labels_2021 == expected_labels_2021
assert {"mol:13CO", "mol:C18O", "mol:DCO+", "mol:DNC"}.issubset(
    linked_labels_2021
)

with TemporaryDirectory() as tmp:
    output_dir = Path(tmp)
    assert write_ppd_table(view_2021, output_dir) == fragments_2021
    assert (output_dir / "ppd_table.tex").read_text() == content_2021

assert ppd_table_atom_groups(view_2026) == [(2, 3, 4, 5, 6), (7, 9, 12)]

detections_2026 = ppd_table_detections(view_2026)
assert len(detections_2026) == counts_2026["detections"]
assert (
    sum(detection.molecule.isotopologue_of is not None for detection in detections_2026)
    == counts_2026["isotopologue_detections"]
)
assert all(detection.status == "secure" for detection in detections_2026)

content_2026 = ppd_table_fragments(view_2026)["ppd_table.tex"]
linked_labels_2026 = re.findall(r"\\molref\{(mol:[^}]+)\}", content_2026)
assert len(linked_labels_2026) == counts_2026["linked_labels"]
assert len(linked_labels_2026) == len(set(linked_labels_2026))
assert set(linked_labels_2026) == {
    detection.molecule.label
    for detection in detections_2026
}
assert {"mol:34SO", "mol:33SO", "mol:H213CO", "mol:c-C2H4O", "mol:HC5N"}.issubset(
    linked_labels_2026
)

print("LaTeX PPD table verification passed")
