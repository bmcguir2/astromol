from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from astromol.census import CensusView
from astromol.database import Database
from astromol.figures import (
    ASTROMOL_BLUE,
    FIGURE_AXES_BOUNDS,
    FIGURE_SIZE,
    NRAO_BLUE,
    detection_rate_by_atoms_data,
    plot_detection_rate_by_atoms,
    plot_detection_rate_by_atoms_comparison,
    write_detection_rate_by_atoms_comparison_plot,
    write_detection_rate_by_atoms_plot,
)


def values_by_label(data):
    """Return count, first year, and rate keyed by point label."""
    return {
        point.label: (point.count, point.first_year, point.rate)
        for point in data.points
    }


db = Database()
view_2021 = CensusView.for_census(db, "2021")
view_2026 = CensusView.for_census(db, "2026")

data_2021 = detection_rate_by_atoms_data(view_2021)
assert data_2021.end_year == 2021
values_2021 = values_by_label(data_2021)
assert values_2021["2"][:2] == (41, 1937)
assert np.isclose(values_2021["2"][2], 41 / (2021 - 1937 + 1))
assert values_2021["3"][:2] == (45, 1969)
assert np.isclose(values_2021["3"][2], 45 / (2021 - 1969 + 1))
assert values_2021["12"][:2] == (5, 2001)
assert np.isclose(values_2021["12"][2], 5 / (2021 - 2001 + 1))
assert values_2021["13+"][:2] == (2, 2018)
assert values_2021["PAHs"][:2] == (3, 2021)
assert values_2021["Fullerenes"][:2] == (3, 2010)

visible_labels = [point.label for point in data_2021.visible_points]
assert visible_labels == [str(natoms) for natoms in range(2, 13)]

data_2026 = detection_rate_by_atoms_data(view_2026)
values_2026 = values_by_label(data_2026)
assert data_2026.end_year == 2026
assert values_2026["13+"][:2] == (7, 2018)
assert np.isclose(values_2026["13+"][2], 7 / (2026 - 2018 + 1))
assert values_2026["PAHs"][:2] == (9, 2021)
assert np.isclose(values_2026["PAHs"][2], 9 / (2026 - 2021 + 1))

with TemporaryDirectory() as tmp:
    output_dir = Path(tmp)
    figure, axes = plot_detection_rate_by_atoms(data_2021)
    assert axes.get_xlim() == (1.0, 12.8)
    assert axes.get_ylim() == (0.0, 1.0)
    assert len(axes.collections) == 1
    assert len(axes.collections[0].get_offsets()) == 11
    assert len(axes.texts) == 2

    output_path = write_detection_rate_by_atoms_plot(
        data_2021,
        output_dir / "rate_by_atoms.pdf",
    )
    assert output_path.exists()
    assert output_path.stat().st_size > 0

    comparison_figure, comparison_axes = plot_detection_rate_by_atoms_comparison(
        data_2026,
        data_2021,
        current_label="2026",
    )
    assert np.allclose(comparison_figure.get_size_inches(), FIGURE_SIZE)
    assert np.allclose(
        comparison_axes.get_position().bounds,
        FIGURE_AXES_BOUNDS,
    )
    assert comparison_axes.get_xlim() == (1.0, 12.8)
    assert comparison_axes.get_ylim() == (0.0, 1.0)
    assert len(comparison_axes.collections) == 2
    assert len(comparison_axes.collections[0].get_offsets()) == 11
    assert len(comparison_axes.collections[1].get_offsets()) == 11
    assert comparison_axes.texts[1].get_text().startswith("Marker area")

    baseline_sizes = comparison_axes.collections[0].get_sizes()
    current_sizes = comparison_axes.collections[1].get_sizes()
    assert np.all(current_sizes >= baseline_sizes)
    assert np.isclose(
        current_sizes[-1] - baseline_sizes[-1],
        (values_2026["12"][0] - values_2021["12"][0]) * 70,
    )
    assert np.allclose(
        comparison_axes.collections[1].get_facecolors()[0, :3],
        [30 / 255, 144 / 255, 255 / 255],
    )
    assert np.allclose(
        comparison_axes.collections[1].get_edgecolors()[0, :3],
        [10 / 255, 21 / 255, 137 / 255],
    )
    assert ASTROMOL_BLUE == "dodgerblue"
    assert NRAO_BLUE == "#0A1589"

    comparison_output_path = write_detection_rate_by_atoms_comparison_plot(
        data_2026,
        data_2021,
        output_dir / "rate_by_atoms_comparison.pdf",
        current_label="2026",
    )
    assert comparison_output_path.exists()
    assert comparison_output_path.stat().st_size > 0

    import matplotlib.pyplot as plt

    plt.close(figure)
    plt.close(comparison_figure)

print("Detection-rate-by-atoms figure verification passed")
