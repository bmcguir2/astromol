"""Independent fixtures for defects found in the repository review."""

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from astromol.census import CensusView
from astromol.database import Database, DATA_DIR
from astromol.models import Detection, Molecule, RecordHistory, DipoleMoment, RotationalConstants
from astromol.slides import SlideMoleculeEntry, selected_records_last_modified
from astromol.validation import validate_database
from baseline import load_production_baseline


def tiny_db(*, dipole=None, rotcon=None, accepted=True):
    mol = Molecule(name="example", formula="CO", label="mol:EXAMPLE", dipole=DipoleMoment(**dipole) if dipole else None, rotcon=RotationalConstants(**rotcon) if rotcon else None)
    det = Detection(id="det:EXAMPLE:ism-csm:2024", molecule=mol, sources=[], telescopes=[],
                    wavelengths=["mm"], year=2024, type="ISM/CSM",
                    history=RecordHistory(introduced={"date": "2024-01-01"},
                                          accepted={"date": "2024-01-01", "context": "confirmed"} if accepted else None))
    return SimpleNamespace(molecules={mol.label: mol}, detections=[det],
                           detections_by_id={det.id: det}, sources={}, telescopes={})


@pytest.mark.parametrize("value", [float("nan"), float("inf"), -float("inf")])
def test_nonfinite_spectroscopy_is_rejected(value):
    report = validate_database(tiny_db(dipole={"a": value}, rotcon={"A": value}))
    assert {issue.code for issue in report.errors} >= {"dipole-nonfinite", "rotcon-nonfinite"}


def test_invalid_rotor_secure_acceptance_and_self_relationship_are_rejected():
    db = tiny_db(rotcon={"A": 1, "B": 5, "C": 2}, accepted=False)
    det = db.detections[0]
    det.confirms = [det.id]
    det.confirmed_by = [det.id]
    assert {issue.code for issue in validate_database(db).errors} >= {
        "rotcon-order", "history-secure-acceptance", "relationship-self"}


def test_unpublished_year_alias_includes_date_only_tentative_records():
    db = tiny_db(accepted=False)
    db.detections[0].status = "tentative"
    assert CensusView.for_census(db, "2026").is_current
    assert CensusView.for_census(db, "2026").detections(include_tentative=True) == db.detections
    assert not CensusView.for_census(db, "2021").detections(include_tentative=True)


def test_production_membership_against_raw_json_and_reviewed_snapshot():
    raw = json.loads((DATA_DIR / "detections.json").read_text())
    db = Database()
    snapshot = load_production_baseline()["membership"]
    for scope in ("2018", "2021", "current"):
        expected = []
        for record in raw:
            event = (record.get("history") or {}).get("accepted" if record.get("status", "secure") == "secure" else "introduced")
            if not event:
                continue
            if scope == "current" or (event.get("census") is not None and int(event["census"]) <= int(scope)):
                expected.append(record["id"])
        view = CensusView.current(db) if scope == "current" else CensusView.for_census(db, scope)
        actual = sorted(d.id for d in view.detections(include_tentative=True, include_disputed=True, include_isotopologues=True))
        assert actual == sorted(expected) == snapshot[scope]


def test_slide_date_includes_detection_and_facility_updates():
    db = tiny_db()
    mol, det = next(iter(db.molecules.values())), db.detections[0]
    mol.history = RecordHistory(last_modified="2026-05-07")
    det.history.last_modified = "2026-07-07"
    entries = [SlideMoleculeEntry(mol, "CO")]
    assert selected_records_last_modified(entries, [det]) == "2026-07-07"
    det.telescopes = [SimpleNamespace(history=RecordHistory(last_modified="2026-08-01"))]
    assert selected_records_last_modified(entries, [det]) == "2026-08-01"


@pytest.mark.parametrize("exact_constant", [False, True])
def test_explicit_endpoint_controls_both_figure_and_scalar_rates(monkeypatch, exact_constant):
    import numpy as np
    from astromol.figures import cumulative_detection_data, linear_cumulative_rate, scopes_by_year_data
    from astromol.latex import detection_rate_since
    if exact_constant:
        # NumPy may trim a zero slope, leaving only the constant coefficient.
        monkeypatch.setattr(np.polynomial.Polynomial, "fit", lambda *args, **kwargs: np.polynomial.Polynomial([1.0]))
    db = tiny_db()
    view = CensusView.current(db, end_year=2027)
    assert cumulative_detection_data(view, start_year=2024).end_year == 2027
    assert detection_rate_since(view, 2024) == pytest.approx(0.0, abs=1e-12)
    assert linear_cumulative_rate(np.arange(2024, 2028), np.ones(4), start_year=2024) == pytest.approx(0.0, abs=1e-12)
    telescope = SimpleNamespace(nick="facility", name="Facility", shortname="F", latex_name=None, built=2024)
    db.telescopes[telescope.nick] = telescope
    db.detections[0].telescopes = [telescope]
    assert scopes_by_year_data(view, start_year=2024, min_detections=1).series[0].rate == pytest.approx(0.0, abs=1e-12)


def test_explicit_figure_exports_exclude_incidental_imports():
    import astromol.figures as figures
    assert "write_cumulative_detections_plot" in figures.__all__
    assert not hasattr(figures, "np")
    assert not hasattr(figures, "dataclass")


def test_shared_counts_preserve_molecule_and_detection_filters():
    from collections import Counter
    from copy import deepcopy
    db = tiny_db()
    first = db.detections[0]
    first.sources = [SimpleNamespace(nick="carbon-source", name="Carbon", latex_name=None, type="Diffuse Cloud")]
    first.telescopes = [SimpleNamespace(nick="carbon-telescope", name="Carbon", shortname="C", latex_name=None)]
    second = deepcopy(first)
    second.id = "det:WATER:ism-csm:2025"
    second.molecule = Molecule(name="water", formula="H2O", label="mol:WATER")
    second.sources[0].nick = "water-source"
    second.telescopes[0].nick = "water-telescope"
    db.detections.append(second)
    view = CensusView.current(db)
    for filtered in (view.filtered(molecule_filter=lambda m: m.label == "mol:EXAMPLE"),
                     view.filtered(detection_filter=lambda d: d.id == first.id)):
        assert filtered.source_counts() == Counter({"carbon-source": 1})
        assert filtered.facility_counts() == Counter({"carbon-telescope": 1})
        assert filtered.source_counts(group_diffuse_cloud=True) == Counter({"DiffuseCloud": 1})
        assert filtered.facility_counts(key="latex_name") == Counter({"C": 1})
