import re
from pathlib import Path
from tempfile import TemporaryDirectory

from astromol.census import CensusView
from astromol.database import Database
from astromol.latex import (
    ice_table_detections,
    ice_table_fragments,
    write_ice_table,
)


db = Database()
view_2021 = CensusView.for_census(db, "2021")
view_2026 = CensusView.for_census(db, "2026")
current_view = CensusView.current(db)

detections_2021 = ice_table_detections(view_2021)
assert len(detections_2021) == 10
assert all(detection.molecule.isotopologue_of is None for detection in detections_2021)
assert sum(detection.status == "tentative" for detection in detections_2021) == 1
assert all(detection.refs.get("observation") for detection in detections_2021)

fragments_2021 = ice_table_fragments(view_2021)
assert set(fragments_2021) == {"ice_table.tex"}
content_2021 = fragments_2021["ice_table.tex"]
assert content_2021.count(r"\molref{mol:") == 10
assert content_2021.count(r"\citet{") == 11
assert content_2021.count(r"$^{\dagger}$") == 2
assert r"\molref{mol:OCN-}{OCN-}$^{\dagger}$" in content_2021

expected_first_row = r"\molref{mol:CO}{CO}" + "\t&\t1" + r"\\"
assert content_2021.splitlines()[8] == expected_first_row

linked_labels_2021 = set(re.findall(r"\\molref\{(mol:[^}]+)\}", content_2021))
expected_labels_2021 = {detection.molecule.label for detection in detections_2021}
assert linked_labels_2021 == expected_labels_2021
assert linked_labels_2021 == {
    "mol:CH3OH",
    "mol:CH4",
    "mol:CO",
    "mol:CO2",
    "mol:H2CO",
    "mol:H2O",
    "mol:HCOOH",
    "mol:NH3",
    "mol:OCN-",
    "mol:OCS",
}

with TemporaryDirectory() as tmp:
    output_dir = Path(tmp)
    assert write_ice_table(view_2021, output_dir) == fragments_2021
    assert (output_dir / "ice_table.tex").read_text() == content_2021

detections_2026 = ice_table_detections(view_2026)
assert len(detections_2026) == 15
assert all(detection.molecule.isotopologue_of is None for detection in detections_2026)
assert sum(detection.status == "tentative" for detection in detections_2026) == 1

content_2026 = ice_table_fragments(view_2026)["ice_table.tex"]
linked_labels_2026 = re.findall(r"\\molref\{(mol:[^}]+)\}", content_2026)
assert len(linked_labels_2026) == 15
assert len(linked_labels_2026) == len(set(linked_labels_2026))
assert set(linked_labels_2026) == {
    detection.molecule.label
    for detection in detections_2026
}
assert "mol:OCN-" in linked_labels_2026
assert r"\molref{mol:OCN-}{OCN-}$^{\dagger}$" in content_2026
assert content_2026.count(r"\citet{") == 12

tentative_current = ice_table_detections(current_view)
assert len(tentative_current) == 15
assert tentative_current[-1].molecule.label == "mol:OCN-"
assert tentative_current[-1].status == "tentative"
tentative_content = ice_table_fragments(current_view)["ice_table.tex"]
assert r"\molref{mol:OCN-}{OCN-}$^{\dagger}$" in tentative_content

secure_only_current = ice_table_detections(current_view, include_tentative=False)
assert len(secure_only_current) == 14
assert all(detection.status == "secure" for detection in secure_only_current)
secure_only_content = ice_table_fragments(
    current_view,
    include_tentative=False,
)["ice_table.tex"]
assert "mol:OCN-" not in secure_only_content
assert r"$^{\dagger}$" not in secure_only_content

print("LaTeX ice table verification passed")
