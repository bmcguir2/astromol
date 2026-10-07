"""Exercise the APIs that establish our advertised dependency floors."""

import matplotlib.pyplot as plt

from astromol.census import CensusView
from astromol.database import Database
from astromol.figures import mass_by_wavelength_data, plot_mass_by_wavelength_boxplot


def test_mass_and_horizontal_boxplot():
    db = Database()
    assert 27 < db.get_molecule("mol:CO").mass < 29
    fig, ax = plot_mass_by_wavelength_boxplot(mass_by_wavelength_data(CensusView.current(db)))
    try:
        assert "Mass" in ax.get_xlabel()
        assert ax.patches or ax.lines
    finally:
        plt.close(fig)
