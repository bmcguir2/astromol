import re
from pathlib import Path
from tempfile import TemporaryDirectory

from astromol.census import CensusView
from astromol.database import Database
from astromol.latex import (
    exoplanet_table_detections,
    exoplanet_table_fragments,
    write_exoplanet_table,
)
from baseline import load_production_baseline


db = Database()
view_2021 = CensusView.for_census(db, "2021")
view_current = CensusView.current(db, end_year=load_production_baseline()["analysis_end_year"])
counts_current = load_production_baseline()["regression_counts"][
    "latex_exoplanet_table_current"
]

detections_2021 = exoplanet_table_detections(view_2021)
assert len(detections_2021) == 9
assert all(detection.molecule.isotopologue_of is None for detection in detections_2021)
assert all(detection.status == "secure" for detection in detections_2021)
assert all(detection.refs.get("observation") for detection in detections_2021)

fragments_2021 = exoplanet_table_fragments(view_2021)
assert set(fragments_2021) == {"exo_table.tex"}
content_2021 = fragments_2021["exo_table.tex"]
assert content_2021.count(r"\molref{mol:") == 9
assert content_2021.count(r"\citet{") == 19
assert r"$^{\dagger}$" not in content_2021

expected_first_row = r"\molref{mol:OH}{OH}" + "\t&\t1" + r"\\"
assert content_2021.splitlines()[8] == expected_first_row

linked_labels_2021 = set(re.findall(r"\\molref\{(mol:[^}]+)\}", content_2021))
expected_labels_2021 = {detection.molecule.label for detection in detections_2021}
assert linked_labels_2021 == expected_labels_2021
assert linked_labels_2021 == {
    "mol:C2H2",
    "mol:CH4",
    "mol:CO",
    "mol:CO2",
    "mol:H2O",
    "mol:HCN",
    "mol:NH3",
    "mol:OH",
    "mol:TiO",
}

with TemporaryDirectory() as tmp:
    output_dir = Path(tmp)
    assert write_exoplanet_table(view_2021, output_dir) == fragments_2021
    assert (output_dir / "exo_table.tex").read_text() == content_2021

detections_current = exoplanet_table_detections(view_current)
assert len(detections_current) == counts_current["detections"]
assert all(detection.molecule.isotopologue_of is None for detection in detections_current)

content_current = exoplanet_table_fragments(view_current)["exo_table.tex"]
linked_labels_current = re.findall(r"\\molref\{(mol:[^}]+)\}", content_current)
assert len(linked_labels_current) == counts_current["linked_labels"]
assert len(linked_labels_current) == len(set(linked_labels_current))
assert set(linked_labels_current) == {
    detection.molecule.label
    for detection in detections_current
}
assert "mol:13CO" not in linked_labels_current
assert "mol:CH3D" not in linked_labels_current

isotope_detections_current = exoplanet_table_detections(
    view_current,
    include_isotopologues=True,
)
assert len(isotope_detections_current) == counts_current["detections_with_isotopologues"]
assert sum(
    detection.molecule.isotopologue_of is not None
    for detection in isotope_detections_current
) == counts_current["isotopologue_detections"]

isotope_content_current = exoplanet_table_fragments(
    view_current,
    include_isotopologues=True,
)["exo_table.tex"]
assert r"\molref{mol:13CO}{^{13}CO}" in isotope_content_current
assert r"\molref{mol:CH3D}{CH3D}" in isotope_content_current

print("LaTeX exoplanet table verification passed")
