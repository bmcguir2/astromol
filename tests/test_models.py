import pytest

from astromol.models import DipoleMoment, Molecule, RotationalConstants


def test_dipole_total_returns_none_for_unknown_placeholder():
    assert DipoleMoment(a="*").total is None


def test_dipole_total_returns_none_for_partially_unknown_components():
    assert DipoleMoment(a=1.6, b="*").total is None


def test_dipole_total_uses_numeric_components():
    assert DipoleMoment(a=3.0, b=4.0).total == 5.0


def test_isotope_formula_fallback_properties():
    molecule = Molecule(
        name="aluminum monofluoride",
        formula="[26Al]F",
        label="mol:26AlF",
    )

    assert molecule.atom_counts == {"Al": 1, "F": 1}
    assert molecule.isotope_counts == {"26Al": 1, "F": 1}
    assert molecule.natoms == 2
    assert molecule.nominal_mass == 45
    assert molecule.mass == pytest.approx(44.98529504273)


def test_charge_and_radical_properties():
    cation = Molecule(name="nitric oxide cation", formula="NO+", label="mol:NO+")

    assert cation.charge == 1
    assert cation.cation
    assert not cation.anion
    assert not cation.neutral
    assert cation.nelectrons == 14
    assert not cation.radical

    radical = Molecule(name="methylidyne", formula="CH", label="mol:CH")
    assert radical.odd_electron
    assert radical.radical

    override = Molecule(
        name="curated non-radical",
        formula="CH",
        label="mol:CH-nonradical",
        radical_override=False,
    )
    assert not override.radical


def test_degree_of_unsaturation_properties():
    acetonitrile = Molecule(
        name="acetonitrile",
        formula="C2H3N",
        label="mol:CH3CN",
    )
    assert acetonitrile.du == 2.0
    assert acetonitrile.maxdu == 3.5

    silicon_dicarbide = Molecule(
        name="silicon dicarbide",
        formula="SiC2",
        label="mol:c-SiC2",
    )
    assert silicon_dicarbide.du is None
    assert silicon_dicarbide.maxdu is None


def test_kappa_properties():
    linear = Molecule(
        name="carbon monoxide",
        formula="CO",
        label="mol:CO",
        rotcon=RotationalConstants(B=57635.968),
    )
    assert linear.is_linear
    assert linear.kappa == -1.0

    asymmetric_top = Molecule(
        name="water",
        formula="H2O",
        label="mol:H2O",
        rotcon=RotationalConstants(A=10.0, B=5.0, C=2.0),
    )
    assert not asymmetric_top.is_linear
    assert asymmetric_top.kappa == pytest.approx(-0.25)
