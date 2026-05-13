from __future__ import annotations

from dataclasses import dataclass

from astromol.census import CensusView
from astromol.models import Detection, Molecule, RecordHistory


@dataclass
class TinyDatabase:
    molecules: dict[str, Molecule]
    detections: list[Detection]


def history(
    *,
    introduced_census: str,
    accepted_census: str | None,
    introduced_context: str = "confirmed",
    accepted_context: str = "confirmed",
) -> RecordHistory:
    accepted = None
    if accepted_census is not None:
        accepted = {
            "census": accepted_census,
            "context": accepted_context,
        }

    return RecordHistory(
        introduced={
            "census": introduced_census,
            "context": introduced_context,
        },
        accepted=accepted,
        last_modified="2026-05-13",
        events=[],
    )


def molecule(
    label: str,
    formula: str,
    *,
    introduced_census: str,
    accepted_census: str | None,
    introduced_context: str = "confirmed",
    isotopologue_of: str | None = None,
) -> Molecule:
    return Molecule(
        name=label,
        formula=formula,
        label=label,
        isotopologue_of=isotopologue_of,
        history=history(
            introduced_census=introduced_census,
            accepted_census=accepted_census,
            introduced_context=introduced_context,
        ),
    )


def detection(
    detection_id: str,
    mol: Molecule,
    *,
    introduced_census: str,
    accepted_census: str | None,
    status: str = "secure",
    detection_type: str = "ISM/CSM",
) -> Detection:
    introduced_context = "confirmed" if status == "secure" else status
    accepted_context = f"confirmed_{detection_type.lower().replace('/', '_')}"

    return Detection(
        id=detection_id,
        molecule=mol,
        sources=[],
        telescopes=[],
        wavelengths=["mm"],
        year=int(introduced_census),
        type=detection_type,
        status=status,
        history=history(
            introduced_census=introduced_census,
            accepted_census=accepted_census,
            introduced_context=introduced_context,
            accepted_context=accepted_context,
        ),
    )


def tiny_database() -> TinyDatabase:
    parent = molecule(
        "mol:CO",
        "CO",
        introduced_census="2018",
        accepted_census="2021",
    )
    isotope = molecule(
        "mol:13CO",
        "[13C]O",
        introduced_census="2021",
        accepted_census="2021",
        isotopologue_of="mol:CO",
    )
    future = molecule(
        "mol:C2H5OH",
        "C2H6O",
        introduced_census="2026",
        accepted_census="2027",
    )
    tentative = molecule(
        "mol:tentative",
        "C2H4",
        introduced_census="2021",
        accepted_census=None,
        introduced_context="tentative",
    )
    disputed = molecule(
        "mol:disputed",
        "C3H",
        introduced_census="2018",
        accepted_census=None,
        introduced_context="disputed",
    )

    detections = [
        detection(
            "det:CO:ism-csm:2021",
            parent,
            introduced_census="2021",
            accepted_census="2021",
        ),
        detection(
            "det:13CO:isotopologue:2021",
            isotope,
            introduced_census="2021",
            accepted_census="2021",
        ),
        detection(
            "det:C2H5OH:ism-csm:2027",
            future,
            introduced_census="2026",
            accepted_census="2027",
        ),
        detection(
            "det:tentative:ism-csm:2021",
            tentative,
            introduced_census="2021",
            accepted_census=None,
            status="tentative",
        ),
        detection(
            "det:disputed:ism-csm:2018",
            disputed,
            introduced_census="2018",
            accepted_census=None,
            status="disputed",
        ),
    ]

    molecules = {
        mol.label: mol
        for mol in [parent, isotope, future, tentative, disputed]
    }
    return TinyDatabase(molecules=molecules, detections=detections)


def detection_ids(detections: list[Detection]) -> set[str]:
    return {detection.id for detection in detections}


def molecule_labels(molecules: list[Molecule]) -> set[str]:
    return {molecule.label for molecule in molecules}


def test_census_boundary_uses_accepted_census_for_secure_records():
    db = tiny_database()

    assert detection_ids(CensusView.for_census(db, "2018").ism_detections()) == set()
    assert detection_ids(CensusView.for_census(db, "2021").ism_detections()) == {
        "det:CO:ism-csm:2021"
    }
    assert detection_ids(CensusView.for_census(db, "2026").ism_detections()) == {
        "det:CO:ism-csm:2021"
    }
    assert detection_ids(CensusView.current(db).ism_detections()) == {
        "det:CO:ism-csm:2021",
        "det:C2H5OH:ism-csm:2027",
    }


def test_tentative_and_disputed_records_use_introduced_census_when_requested():
    db = tiny_database()
    view_2021 = CensusView.for_census(db, "2021")

    assert "det:tentative:ism-csm:2021" not in detection_ids(view_2021.ism_detections())
    assert "det:disputed:ism-csm:2018" not in detection_ids(view_2021.ism_detections())

    assert "det:tentative:ism-csm:2021" in detection_ids(
        view_2021.ism_detections(include_tentative=True)
    )
    assert "det:disputed:ism-csm:2018" in detection_ids(
        view_2021.ism_detections(include_disputed=True)
    )

    assert "det:tentative:ism-csm:2021" not in detection_ids(
        CensusView.for_census(db, "2018").ism_detections(include_tentative=True)
    )
    assert "det:disputed:ism-csm:2018" in detection_ids(
        CensusView.for_census(db, "2018").ism_detections(include_disputed=True)
    )


def test_isotopologues_are_excluded_by_default_and_included_on_request():
    db = tiny_database()
    view = CensusView.for_census(db, "2021")

    assert detection_ids(view.detections(include_isotopologues=False)) == {
        "det:CO:ism-csm:2021",
    }
    assert detection_ids(view.detections(include_isotopologues=True)) == {
        "det:CO:ism-csm:2021",
        "det:13CO:isotopologue:2021",
    }

    assert molecule_labels(view.accepted_molecules()) == {"mol:CO"}
    assert molecule_labels(view.accepted_molecules(include_isotopologues=True)) == {
        "mol:CO",
        "mol:13CO",
    }


def test_context_molecules_follow_detection_scope_not_molecule_acceptance_only():
    db = tiny_database()
    view = CensusView.for_census(db, "2021")

    assert "mol:tentative" not in molecule_labels(view.ism_molecules())
    assert "mol:tentative" in molecule_labels(
        view.ism_molecules(include_tentative=True)
    )
    assert "mol:tentative" not in molecule_labels(view.accepted_molecules())

