from pathlib import Path
from tempfile import TemporaryDirectory

from astromol.census import CensusView
from astromol.database import Database
from astromol.figures import (
    molecules_by_wavelength_atoms_data,
    plot_molecules_by_wavelength_atoms,
    plot_molecules_by_wavelength_atoms_bubble_heatmap,
    write_molecules_by_wavelength_atoms_bubble_heatmap,
    write_molecules_by_wavelength_atoms_plot,
)


db = Database()
view_2021 = CensusView.for_census(db, "2021")

data_2021 = molecules_by_wavelength_atoms_data(view_2021)
assert data_2021.molecule_count == 237
assert data_2021.counts == {
    "cm": 89,
    "mm": 137,
    "sub-mm": 13,
    "IR": 15,
    "Vis": 2,
    "UV": 6,
}
assert data_2021.max_atoms == 19
matrix_2021 = data_2021.matrix()
assert matrix_2021.shape == (6, 12)
assert matrix_2021[0, 3] == 20  # cm, 5 atoms
assert matrix_2021[1, 1] == 36  # mm, 3 atoms
assert matrix_2021[5, 0] == 6   # UV, 2 atoms

with TemporaryDirectory() as tmp:
    output_dir = Path(tmp)
    figure, axes = plot_molecules_by_wavelength_atoms(data_2021)
    assert len(axes) == 6
    assert axes[0].get_ylabel() == "Probability Density Estimate"
    assert axes[3].get_xlabel() == "# of Atoms"
    assert axes[4].get_ylabel() == "# of Detected Molecules"
    assert [text.get_text() for text in axes[0].texts] == ["cm"]
    assert [text.get_text() for text in axes[5].texts] == ["UV"]

    output_path = write_molecules_by_wavelength_atoms_plot(
        data_2021,
        output_dir / "mols_waves_by_atoms.pdf",
    )
    assert output_path.exists()
    assert output_path.stat().st_size > 0

    bubble_figure, bubble_ax = plot_molecules_by_wavelength_atoms_bubble_heatmap(
        data_2021,
    )
    assert len(bubble_figure.axes) == 2
    assert bubble_ax.get_title() == ""
    assert bubble_ax.get_xlabel() == "# of Atoms"
    assert [label.get_text() for label in bubble_ax.get_yticklabels()] == [
        "cm",
        "mm",
        "sub-mm",
        "IR",
        "Vis",
        "UV",
    ]
    assert any(text.get_text() == "36" for text in bubble_ax.texts)

    bubble_output_path = write_molecules_by_wavelength_atoms_bubble_heatmap(
        data_2021,
        output_dir / "mols_waves_by_atoms_bubble_heatmap.pdf",
    )
    assert bubble_output_path.exists()
    assert bubble_output_path.stat().st_size > 0

    percent_bubble_figure, percent_bubble_ax = (
        plot_molecules_by_wavelength_atoms_bubble_heatmap(
            data_2021,
            value_mode="row_percent",
        )
    )
    assert len(percent_bubble_figure.axes) == 2
    assert any(text.get_text() == "100" for text in percent_bubble_ax.texts)

    percent_bubble_output_path = write_molecules_by_wavelength_atoms_bubble_heatmap(
        data_2021,
        output_dir / "mols_waves_by_atoms_bubble_heatmap_row_percent.pdf",
        value_mode="row_percent",
    )
    assert percent_bubble_output_path.exists()
    assert percent_bubble_output_path.stat().st_size > 0

    import matplotlib.pyplot as plt

    plt.close(figure)
    plt.close(bubble_figure)
    plt.close(percent_bubble_figure)

print("Molecules-by-wavelength atom-count figure verification completed")
