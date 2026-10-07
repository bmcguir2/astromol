import re
from pathlib import Path
from tempfile import TemporaryDirectory

from baseline import load_production_baseline

from astromol.census import CensusView
from astromol.database import Database
from astromol.latex import (
    exgal_table_atom_groups,
    exgal_table_columns,
    exgal_table_detections,
    exgal_table_fragments,
    write_exgal_table,
)


db = Database()
view_2021 = CensusView.for_census(db, "2021")
view_current = CensusView.current(db, end_year=load_production_baseline()["analysis_end_year"])
counts_current = load_production_baseline()["regression_counts"][
    "latex_exgal_table_current"
]

assert exgal_table_atom_groups(view_2021) == [(2, 3, 4, 5), (6, 7, 8, 9, 12)]

columns_2021 = exgal_table_columns(view_2021)
assert [[len(column) for column in group] for group in columns_2021] == [
    [22, 18, 12, 9],
    [5, 5, 2, 1, 1],
]

detections_2021 = exgal_table_detections(view_2021)
assert len(detections_2021) == 75
assert sum(detection.status == "tentative" for detection in detections_2021) == 2
assert all(detection.molecule.isotopologue_of is None for detection in detections_2021)
assert all(detection.refs.get("observation") for detection in detections_2021)

fragments_2021 = exgal_table_fragments(view_2021)
assert set(fragments_2021) == {"exgal_table.tex"}
content_2021 = fragments_2021["exgal_table.tex"]
assert content_2021.count(r"\molref{mol:") == 75
assert content_2021.count(r"\citet{") == 52
assert content_2021.count(r"$^{\dagger}$") == 3
assert r"\molref{mol:c-C3H}{c-C3H}$^{\dagger}$" in content_2021
assert r"\molref{mol:HC5N}{HC5N}$^{\dagger}$" in content_2021

expected_first_row = (
    r"\molref{mol:CH}{CH}"
    "\t&\t1\t&\t"
    r"\molref{mol:H2O}{H2O}"
    "\t&\t2\t&\t"
    r"\molref{mol:NH3}{NH3}"
    "\t&\t3\t&\t"
    r"\molref{mol:HC3N}{HC3N}"
    "\t&\t4, 5"
    r"\\"
)
assert content_2021.splitlines()[9] == expected_first_row

linked_labels_2021 = set(re.findall(r"\\molref\{(mol:[^}]+)\}", content_2021))
expected_labels_2021 = {detection.molecule.label for detection in detections_2021}
assert linked_labels_2021 == expected_labels_2021

with TemporaryDirectory() as tmp:
    output_dir = Path(tmp)
    assert write_exgal_table(view_2021, output_dir) == fragments_2021
    assert (output_dir / "exgal_table.tex").read_text() == content_2021

detections_current = exgal_table_detections(view_current)
assert len(detections_current) == counts_current["detections"]
assert sum(
    detection.status == "secure" for detection in detections_current
) == counts_current["secure_detections"]
assert sum(
    detection.status == "tentative" for detection in detections_current
) == counts_current["tentative_detections"]
assert all(detection.molecule.isotopologue_of is None for detection in detections_current)
assert all(detection.refs.get("observation") for detection in detections_current)

content_current = exgal_table_fragments(view_current)["exgal_table.tex"]
linked_labels_current = re.findall(r"\\molref\{(mol:[^}]+)\}", content_current)
assert len(linked_labels_current) == counts_current["linked_labels"]
assert len(linked_labels_current) == len(set(linked_labels_current))
assert set(linked_labels_current) == {
    detection.molecule.label
    for detection in detections_current
}
observation_bibcodes_current = {
    ref.bibcode
    for detection in detections_current
    for ref in detection.refs.get("observation", [])
}
assert len(observation_bibcodes_current) == counts_current["observation_references"]
assert content_current.count(r"\citet{") == counts_current["observation_references"]
assert content_current.count(r"$^{\dagger}$") == (
    counts_current["tentative_detections"] + 1
)

print("LaTeX exgal table verification passed")
