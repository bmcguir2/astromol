from astromol.models import DipoleMoment


def test_dipole_total_returns_none_for_unknown_placeholder():
    assert DipoleMoment(a="*").total is None


def test_dipole_total_returns_none_for_partially_unknown_components():
    assert DipoleMoment(a=1.6, b="*").total is None


def test_dipole_total_uses_numeric_components():
    assert DipoleMoment(a=3.0, b=4.0).total == 5.0
