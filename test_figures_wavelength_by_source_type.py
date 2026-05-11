from pathlib import Path
from tempfile import TemporaryDirectory

from astromol.census import CensusView
from astromol.database import Database
from astromol.figures import (
    plot_wavelength_by_source_type,
    plot_wavelength_by_source_type_stacked_bar,
    wavelength_by_source_type_data,
    write_wavelength_by_source_type,
    write_wavelength_by_source_type_stacked_bar,
)


db = Database()
view_2021 = CensusView.for_census(db, "2021")
data_2021 = wavelength_by_source_type_data(view_2021)

assert data_2021.molecule_count == 217
assert data_2021.counts == {
    "carbon_star": {"cm": 11, "mm": 46, "sub-mm": 1, "IR": 8, "UV": 0, "Vis": 0},
    "dark_cloud": {"cm": 53, "mm": 25, "sub-mm": 0, "IR": 0, "UV": 0, "Vis": 0},
    "diffuse_cloud": {"cm": 4, "mm": 2, "sub-mm": 6, "IR": 6, "UV": 6, "Vis": 2},
    "sfr": {"cm": 31, "mm": 60, "sub-mm": 5, "IR": 0, "UV": 0, "Vis": 0},
}

view_2026 = CensusView.for_census(db, "2026")
data_2026 = wavelength_by_source_type_data(view_2026)
assert data_2026.molecule_count == 284
assert data_2026.counts == {
    "carbon_star": {"cm": 22, "mm": 51, "sub-mm": 3, "IR": 8, "UV": 0, "Vis": 0},
    "dark_cloud": {"cm": 102, "mm": 29, "sub-mm": 0, "IR": 0, "UV": 0, "Vis": 0},
    "diffuse_cloud": {"cm": 4, "mm": 2, "sub-mm": 6, "IR": 6, "UV": 6, "Vis": 2},
    "sfr": {"cm": 31, "mm": 64, "sub-mm": 6, "IR": 0, "UV": 0, "Vis": 0},
}

with TemporaryDirectory() as tmp:
    output_dir = Path(tmp)

    figure, axes = plot_wavelength_by_source_type(
        data_2021,
        source_label_overrides={"diffuse_cloud": "LOS Clouds"},
    )
    panel_texts = [
        text.get_text()
        for ax in axes
        for text in ax.texts
        if text.get_text()
    ]
    assert "Carbon Stars" in panel_texts
    assert "Dark Clouds" in panel_texts
    assert "LOS Clouds" in panel_texts
    assert "SFRs" in panel_texts
    assert "69.7%" in panel_texts
    assert "67.9%" in panel_texts
    assert "30.8%" in panel_texts
    assert "62.5%" in panel_texts

    legend = axes[1].get_legend()
    assert legend is not None
    assert legend.get_title().get_text() == "Wavelengths"
    assert [
        text.get_text()
        for text in legend.get_texts()
    ] == ["cm", "mm", "sub-mm", "IR", "UV/Vis"]

    output_path = write_wavelength_by_source_type(
        data_2021,
        output_dir / "waves_by_source_type.pdf",
        source_label_overrides={"diffuse_cloud": "LOS Clouds"},
    )
    assert output_path.exists()
    assert output_path.stat().st_size > 0

    stacked_figure, stacked_ax = plot_wavelength_by_source_type_stacked_bar(
        data_2026,
    )
    assert stacked_ax.get_xlabel() == "First-detection wavelength credits (%)"
    assert [
        label.get_text()
        for label in stacked_ax.get_yticklabels()
    ] == [
        "Dark Cloud",
        "Carbon Star",
        "SFR",
        "Diffuse Cloud",
    ]
    assert [
        text.get_text()
        for text in stacked_ax.texts
        if text.get_text().startswith("n=")
    ] == [
        "n=131",
        "n=84",
        "n=101",
        "n=26",
    ]
    stacked_output_path = write_wavelength_by_source_type_stacked_bar(
        data_2026,
        output_dir / "waves_by_source_type_stacked_bar.pdf",
    )
    assert stacked_output_path.exists()
    assert stacked_output_path.stat().st_size > 0

    import matplotlib.pyplot as plt

    plt.close(figure)
    plt.close(stacked_figure)

print("Wavelength-by-source-type figure verification completed")
