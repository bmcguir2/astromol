from astromol.census import CensusView
from astromol.database import Database


def labels(records):
    return {record.molecule.label for record in records}


db = Database()

view_2018 = CensusView.for_census(db, "2018")
assert len(view_2018.ism_molecules()) == 204
assert len(view_2018.exgal_molecules()) == 63
assert len(view_2018.exgal_molecules(include_tentative=True)) == 65
assert len(view_2018.exoplanet_molecules()) == 5
assert len(view_2018.ppd_molecules()) == 23
assert len(view_2018.ppd_molecules(include_isotopologues=True)) == 35

view_2021 = CensusView.for_census(db, "2021")
assert len(view_2021.ism_molecules()) == 240
assert len(view_2021.exgal_molecules()) == 73
assert len(view_2021.exgal_molecules(include_tentative=True)) == 75
assert len(view_2021.exoplanet_molecules()) == 9
assert len(view_2021.ppd_molecules()) == 25
assert len(view_2021.ppd_molecules(include_isotopologues=True)) == 40

source_counts = view_2021.source_counts(
    key="latex_name",
    group_diffuse_cloud=True,
    diffuse_cloud_label="LOS Cloud",
)
assert len(source_counts) == 38
assert sum(source_counts.values()) == 332
assert source_counts["Sgr B2"] == 68
assert source_counts["LOS Cloud"] == 42

facility_counts = view_2021.facility_counts(key="latex_name")
assert len(facility_counts) == 46
assert sum(facility_counts.values()) == 304
assert facility_counts["IRAM 30-m"] == 64
assert facility_counts["NRAO 36-ft"] == 33

view_2026 = CensusView.for_census(db, "2026")
assert len(view_2026.ism_molecules()) == 325
assert len(view_2026.ism_molecules(include_isotopologues=True)) == 335
assert len(view_2026.exoplanet_molecules()) == 11
assert len(view_2026.exoplanet_molecules(include_isotopologues=True)) == 13
assert len(view_2026.ppd_molecules()) == 34
assert len(view_2026.ppd_molecules(include_isotopologues=True)) == 57
current = CensusView.current(db)
for detection_type in ["ISM/CSM", "exgal", "exo", "ppd", "ice"]:
    assert {
        detection.id for detection in view_2026.context_detections(detection_type)
    } == {
        detection.id for detection in current.context_detections(detection_type)
    }
    assert {
        detection.id
        for detection in view_2026.context_detections(
            detection_type,
            include_isotopologues=True,
        )
    } == {
        detection.id
        for detection in current.context_detections(
            detection_type,
            include_isotopologues=True,
        )
    }

print("CensusView verification passed")
