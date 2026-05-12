import os
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from astromol.census import CensusView
from astromol.database import Database
from astromol.figures import (
    COLORBLIND_CUMULATIVE_BY_ATOMS_SERIES,
    FIGURE_AXES_BOUNDS,
    FIGURE_SIZE,
    LEGACY_CUMULATIVE_BY_ATOMS_SERIES,
    ROLLING_RATE_HEATMAP_AXES_BOUNDS,
    ROLLING_RATE_HEATMAP_COLORBAR_BOUNDS,
    ROLLING_RATE_HEATMAP_SIZE,
    cumulative_by_atoms_data,
    plot_cumulative_by_atoms,
    plot_rolling_rate_by_atoms_heatmap,
    plot_stacked_cumulative_by_atoms,
    rolling_rate_by_atoms_heatmap_data,
    write_cumulative_by_atoms_plot,
    write_rolling_rate_by_atoms_heatmap,
    write_stacked_cumulative_by_atoms_plot,
)


def final_counts_by_label(data):
    """Return final cumulative counts keyed by series label."""
    return {series.label: series.final_count for series in data.series}


def colors_by_label(data):
    """Return plotted colors keyed by series label."""
    return {series.label: series.color for series in data.series}


db = Database()
view_2021 = CensusView.for_census(db, "2021")
view_2026 = CensusView.for_census(db, "2026")
view_current = CensusView.current(db)

data_2021 = cumulative_by_atoms_data(view_2021)
assert data_2021.start_year == 1937
assert data_2021.end_year == 2021
assert data_2021.total == 240
assert final_counts_by_label(data_2021) == {
    "2 atoms": 41,
    "3 atoms": 45,
    "4 atoms": 31,
    "5 atoms": 31,
    "6 atoms": 23,
    "7 atoms": 15,
    "8 atoms": 15,
    "9 atoms": 14,
    "10 atoms": 6,
    "11 atoms": 6,
    "12 atoms": 5,
    "13+ atoms": 2,
    "Fullerenes": 3,
    "PAHs": 3,
}
assert colors_by_label(data_2021) == {
    spec["label"]: spec["color"]
    for spec in LEGACY_CUMULATIVE_BY_ATOMS_SERIES
}

data_2026 = cumulative_by_atoms_data(view_2026)
assert data_2026.start_year == 1937
assert data_2026.end_year == 2026
assert data_2026.total == 325
assert final_counts_by_label(data_2026)["13+ atoms"] == 7
assert "mol:c-C6H5CCH" in next(
    series.first_detection_years
    for series in data_2026.series
    if series.label == "13+ atoms"
)
assert final_counts_by_label(data_2026)["PAHs"] == 9
assert final_counts_by_label(data_2026)["Fullerenes"] == 3
assert colors_by_label(data_2026) == {
    spec["label"]: spec["color"]
    for spec in COLORBLIND_CUMULATIVE_BY_ATOMS_SERIES
}

data_current = cumulative_by_atoms_data(view_current)
assert colors_by_label(data_current) == colors_by_label(data_2026)

heatmap_2026 = rolling_rate_by_atoms_heatmap_data(view_2026)
assert heatmap_2026.labels == (
    "2",
    "3",
    "4",
    "5",
    "6",
    "7",
    "8",
    "9",
    "10",
    "11",
    "12",
    "13+",
    "Fuller",
    "PAH",
)
assert heatmap_2026.matrix.shape == (14, len(heatmap_2026.years))
assert heatmap_2026.window == 10
assert heatmap_2026.vmax == 2.0
assert np.allclose(
    heatmap_2026.matrix[:, -1],
    [
        0.8,
        0.9,
        1.1,
        1.8,
        1.6,
        1.3,
        1.2,
        0.8,
        1.0,
        0.6,
        0.6,
        0.7,
        0.0,
        0.9,
    ],
)

with TemporaryDirectory() as tmp:
    output_dir = Path(tmp)
    os.environ.setdefault("MPLCONFIGDIR", str(output_dir / "mplconfig"))
    os.environ.setdefault("XDG_CACHE_HOME", str(output_dir / "cache"))

    figure, axes = plot_cumulative_by_atoms(data_2021)
    assert axes.get_ylim() == (-1.0, 47.25)
    assert len(axes.lines) == 14
    assert len(axes.texts) == 14

    stacked_figure, stacked_axes = plot_stacked_cumulative_by_atoms(data_2021)
    assert np.allclose(stacked_axes.get_ylim(), (0.0, 249.6))
    assert len(stacked_axes.collections) == 14
    assert len(stacked_axes.texts) == 14
    assert np.allclose(figure.get_size_inches(), FIGURE_SIZE)
    assert np.allclose(stacked_figure.get_size_inches(), FIGURE_SIZE)
    assert np.allclose(axes.get_position().bounds, FIGURE_AXES_BOUNDS)
    assert np.allclose(stacked_axes.get_position().bounds, FIGURE_AXES_BOUNDS)
    assert axes.get_xlim() == stacked_axes.get_xlim()
    assert axes.texts[0].get_position() == stacked_axes.texts[0].get_position()
    assert axes.texts[-1].get_position() == stacked_axes.texts[-1].get_position()

    heatmap_figure, heatmap_axes = plot_rolling_rate_by_atoms_heatmap(
        heatmap_2026,
    )
    assert np.allclose(heatmap_figure.get_size_inches(), ROLLING_RATE_HEATMAP_SIZE)
    assert np.allclose(
        heatmap_axes.get_position().bounds,
        ROLLING_RATE_HEATMAP_AXES_BOUNDS,
    )
    assert np.allclose(
        heatmap_figure.axes[1].get_position().bounds,
        ROLLING_RATE_HEATMAP_COLORBAR_BOUNDS,
    )
    assert len(heatmap_axes.images) == 1
    assert len(heatmap_axes.lines) == 0
    assert [label.get_text() for label in heatmap_axes.get_yticklabels()] == list(
        heatmap_2026.labels
    )

    output_path = write_cumulative_by_atoms_plot(
        data_2021,
        output_dir / "cumulative_by_atoms.pdf",
    )
    assert output_path.exists()
    assert output_path.stat().st_size > 0

    stacked_output_path = write_stacked_cumulative_by_atoms_plot(
        data_2021,
        output_dir / "stacked_cumulative_by_atoms.pdf",
    )
    assert stacked_output_path.exists()
    assert stacked_output_path.stat().st_size > 0

    heatmap_output_path = write_rolling_rate_by_atoms_heatmap(
        heatmap_2026,
        output_dir / "rolling_rate_heatmap.pdf",
    )
    assert heatmap_output_path.exists()
    assert heatmap_output_path.stat().st_size > 0

    import matplotlib.pyplot as plt

    plt.close(figure)
    plt.close(stacked_figure)
    plt.close(heatmap_figure)

print("Cumulative detections by atom count figure verification passed")
