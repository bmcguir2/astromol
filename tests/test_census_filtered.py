from astromol.census import CensusView
from astromol.database import Database


def test_filtered_view_restricts_molecules_and_counts():
    db = Database()
    view = CensusView.for_census(db, "2026")
    carbon_view = view.filtered(
        molecule_filter=lambda molecule: molecule.atom_counts.get("C", 0) > 0
    )

    molecules = carbon_view.ism_molecules()

    assert molecules
    assert len(molecules) < len(view.ism_molecules())
    assert all(molecule.atom_counts.get("C", 0) > 0 for molecule in molecules)
    assert sum(carbon_view.source_counts().values()) <= sum(
        view.source_counts().values()
    )


def test_filtered_view_restricts_detections():
    db = Database()
    view = CensusView.for_census(db, "2026")
    mm_view = view.filtered(
        detection_filter=lambda detection: "mm" in detection.wavelengths
    )

    detections = mm_view.ism_detections()

    assert detections
    assert len(detections) < len(view.ism_detections())
    assert all("mm" in detection.wavelengths for detection in detections)
