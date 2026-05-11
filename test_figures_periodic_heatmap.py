from pathlib import Path
from tempfile import TemporaryDirectory

from astromol.census import CensusView
from astromol.database import Database
from astromol.figures import (
    periodic_heatmap_data,
    plot_periodic_heatmap,
    write_periodic_heatmap,
)


db = Database()
view_2021 = CensusView.for_census(db, "2021")

data_2021 = periodic_heatmap_data(view_2021)
assert data_2021.molecule_count == 240
assert data_2021.detected_element_count == 19
assert data_2021.element_counts["C"] == 188
assert data_2021.element_counts["H"] == 173
assert data_2021.element_counts["N"] == 92
assert data_2021.element_counts["O"] == 75
assert data_2021.element_counts["S"] == 30
assert data_2021.element_counts["Si"] == 13

cell_by_symbol = {cell.symbol: cell for cell in data_2021.cells}
assert cell_by_symbol["C"].count == 188
assert cell_by_symbol["C"].group == 14
assert cell_by_symbol["C"].y_position == 5.75
assert cell_by_symbol["He"].count == 1
assert cell_by_symbol["Ne"].count == 0

with TemporaryDirectory() as tmp:
    output_dir = Path(tmp)
    figure, ax = plot_periodic_heatmap(data_2021)
    assert ax.get_xlim() == (0.0, 18.0)
    assert ax.get_ylim() == (0.0, 8.0)
    assert len(ax.patches) == len(data_2021.cells)
    assert any(text.get_text() == "188" for text in ax.texts)

    output_path = write_periodic_heatmap(
        data_2021,
        output_dir / "periodic_heatmap.pdf",
    )
    assert output_path.exists()
    assert output_path.stat().st_size > 0

    import matplotlib.pyplot as plt

    plt.close(figure)

print("Periodic heatmap figure verification passed")
