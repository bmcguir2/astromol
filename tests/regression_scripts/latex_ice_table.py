import re
from pathlib import Path
from tempfile import TemporaryDirectory

from baseline import load_production_baseline

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
counts_2026 = load_production_baseline()["regression_counts"][
    "latex_ice_table_2026"
]


def observation_bibcodes(detections):
    return {
        ref.bibcode
        for detection in detections
        for ref in detection.refs.get("observation", [])
    }


def assert_current_table_invariants(content, detections, expected_counts):
    linked_labels = re.findall(r"\\molref\{(mol:[^}]+)\}", content)
    assert len(detections) == expected_counts["detections"]
    assert len(linked_labels) == expected_counts["linked_labels"]
    assert len(linked_labels) == len(set(linked_labels))
    assert set(linked_labels) == {
        detection.molecule.label
        for detection in detections
    }

    secure_count = sum(detection.status == "secure" for detection in detections)
    tentative_count = sum(
        detection.status == "tentative" for detection in detections
    )
    assert secure_count == expected_counts["secure_detections"]
    assert tentative_count == expected_counts["tentative_detections"]
    assert all(detection.refs.get("observation") for detection in detections)
    assert len(observation_bibcodes(detections)) == expected_counts[
        "observation_references"
    ]
    assert content.count(r"\citet{") == expected_counts["observation_references"]

    expected_dagger_markers = tentative_count + bool(tentative_count)
    assert content.count(r"$^{\dagger}$") == expected_dagger_markers

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
assert all(detection.molecule.isotopologue_of is None for detection in detections_2026)

content_2026 = ice_table_fragments(view_2026)["ice_table.tex"]
assert_current_table_invariants(content_2026, detections_2026, counts_2026)

ocn_2026 = next(
    detection
    for detection in detections_2026
    if detection.molecule.label == "mol:OCN-"
)
assert ocn_2026.id == "det:OCN-:ice:2024"
assert ocn_2026.status == "secure"
assert ocn_2026.confirms == ["det:OCN-:ice:2005"]
assert r"\molref{mol:OCN-}{OCN-}" in content_2026
assert r"\molref{mol:OCN-}{OCN-}$^{\dagger}$" not in content_2026

detections_current = ice_table_detections(current_view)
assert [detection.id for detection in detections_current] == [
    detection.id for detection in detections_2026
]
current_content = ice_table_fragments(current_view)["ice_table.tex"]
assert_current_table_invariants(current_content, detections_current, counts_2026)

secure_only_current = ice_table_detections(current_view, include_tentative=False)
assert len(secure_only_current) == counts_2026["secure_only_detections"]
assert all(detection.status == "secure" for detection in secure_only_current)
secure_only_content = ice_table_fragments(
    current_view,
    include_tentative=False,
)["ice_table.tex"]
assert "mol:OCN-" in secure_only_content
assert r"$^{\dagger}$" not in secure_only_content
if counts_2026["tentative_detections"] == 0:
    assert secure_only_current == detections_current
    assert secure_only_content == current_content

print("LaTeX ice table verification passed")
