from pathlib import Path
from tempfile import TemporaryDirectory

from astromol.census import CensusView
from astromol.database import Database
from astromol.figures import (
    FIGURE_TEXT_SIZE,
    mass_by_wavelength_data,
    plot_mass_by_wavelength,
    plot_mass_by_wavelength_boxplot,
    write_mass_by_wavelength_boxplot,
    write_mass_by_wavelength_plot,
)


db = Database()
view_2021 = CensusView.for_census(db, "2021")

data_2021 = mass_by_wavelength_data(view_2021)
assert data_2021.molecule_count == 240
assert data_2021.counts["cm"] == 89
assert data_2021.counts["mm"] == 137
assert data_2021.counts["IR"] == 18

# These are the current curated database values. UV-Vis differs from the
# published 2021 figure label because the legacy plotting function deduplicated
# UV and visible detections by molecular mass, dropping CH+ after CH.
assert data_2021.counts["sub-mm"] == 13
assert data_2021.counts["UV-Vis"] == 6

data_2021_without_fullerenes = mass_by_wavelength_data(
    view_2021,
    include_fullerenes=False,
)
assert data_2021_without_fullerenes.counts["IR"] == 15

with TemporaryDirectory() as tmp:
    output_dir = Path(tmp)
    figure, ax = plot_mass_by_wavelength(data_2021)
    assert ax.get_title() == ""
    assert ax.get_xlabel() == "Atomic Mass (amu)"
    assert ax.get_ylabel() == "Probability Density Estimate"
    count_labels = {
        text.get_text(): text
        for text in ax.texts
        if text.get_text() in {"89", "137", "13", "18", "6"}
    }
    assert count_labels["89"].get_fontsize() == FIGURE_TEXT_SIZE
    assert count_labels["137"].get_fontsize() == FIGURE_TEXT_SIZE
    assert any(
        text.get_text() == "Detection Wavelengths"
        and text.get_fontsize() == FIGURE_TEXT_SIZE
        for text in ax.texts
    )
    assert any(line.get_transform() == ax.transAxes for line in ax.lines)

    auto_figure, auto_ax = plot_mass_by_wavelength(data_2021, label_mode="auto")
    auto_count_labels = {
        text.get_text(): text
        for text in auto_ax.texts
        if text.get_text() in {"89", "137", "13", "18", "6"}
    }
    assert set(auto_count_labels) == {"89", "137", "13", "18", "6"}
    assert auto_count_labels["13"].get_position() != count_labels["13"].get_position()

    output_path = write_mass_by_wavelength_plot(
        data_2021,
        output_dir / "mass_by_wavelengths_kde.pdf",
    )
    assert output_path.exists()
    assert output_path.stat().st_size > 0

    box_figure, box_ax = plot_mass_by_wavelength_boxplot(data_2021_without_fullerenes)
    assert box_ax.get_title() == ""
    assert box_ax.get_xlabel() == "Atomic Mass (amu)"
    assert [label.get_text() for label in box_ax.get_yticklabels()] == [
        "cm",
        "mm",
        "sub-mm",
        "IR",
        "UV-Vis",
    ]
    marker_lines = [
        line for line in box_ax.lines
        if len(line.get_xdata()) == 2 and set(line.get_xdata()) == {80}
    ]
    assert marker_lines

    box_output_path = write_mass_by_wavelength_boxplot(
        data_2021_without_fullerenes,
        output_dir / "mass_by_wavelengths_boxplot.pdf",
    )
    assert box_output_path.exists()
    assert box_output_path.stat().st_size > 0

    import matplotlib.pyplot as plt

    plt.close(figure)
    plt.close(auto_figure)
    plt.close(box_figure)

print("Mass-by-wavelength figure verification completed with noted mismatches")
