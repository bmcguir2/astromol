from pathlib import Path
from tempfile import TemporaryDirectory

from astromol.census import CensusView
from astromol.database import Database
from astromol.figures import (
    du_histogram_data,
    plot_du_bar_chart,
    plot_du_histogram,
    write_du_bar_chart,
    write_du_histogram,
)


db = Database()
view_2021 = CensusView.for_census(db, "2021")

data_2021 = du_histogram_data(view_2021)
assert data_2021.molecule_count == 196
assert data_2021.saturated_count == 20
assert data_2021.unsaturated_count == 172
assert data_2021.max_du == 12.0
assert data_2021.histogram_counts().tolist() == [
    20,
    10,
    31,
    18,
    30,
    9,
    15,
    8,
    15,
    4,
    7,
    2,
    10,
    4,
    0,
    2,
    1,
    2,
    2,
    0,
    1,
    0,
    0,
    0,
    1,
]
assert data_2021.counts_by_value[-0.5] == 4
assert data_2021.counts_by_value[1.0] == 31

view_2026 = CensusView.for_census(db, "2026")
data_2026 = du_histogram_data(view_2026)
assert data_2026.max_du == 14.0
du_14_labels = [
    label
    for label, value in zip(data_2026.molecule_labels, data_2026.values)
    if value == 14.0
]
assert du_14_labels == [
    "mol:1-C16H9CN",
    "mol:2-C16H9CN",
    "mol:4-C16H9CN",
]

with TemporaryDirectory() as tmp:
    output_dir = Path(tmp)
    figure, ax = plot_du_histogram(data_2021)
    assert ax.get_title() == ""
    assert ax.get_xlabel() == "Degree of Unsaturation"
    assert ax.get_ylabel() == "# of Detected Molecules"
    assert ax.get_ylim() == (0.0, 35.0)
    assert [text.get_text() for text in ax.texts] == [
        "CH$_4$, CH$_3$OH, CH$_3$Cl, ...",
        "HC$_{11}$N",
    ]

    figure_2026, ax_2026 = plot_du_histogram(data_2026)
    assert any(
        text.get_text() == "1-/2-/4-C$_{16}$H$_{9}$CN"
        for text in ax_2026.texts
    )

    bar_figure, bar_ax = plot_du_bar_chart(data_2026)
    assert bar_ax.get_title() == ""
    assert bar_ax.get_xlabel() == "Degree of Unsaturation"
    assert bar_ax.get_ylabel() == "# of Detected Molecules"
    assert not any(text.get_text().isdigit() for text in bar_ax.texts)
    assert [text.get_text() for text in bar_ax.texts] == [
        "CH$_4$, CH$_3$OH, CH$_3$Cl, ...",
        "HC$_{11}$N",
        "1-/2-/4-C$_{16}$H$_{9}$CN",
    ]
    assert min(patch.get_x() for patch in bar_ax.patches) >= -0.2

    output_path = write_du_histogram(
        data_2021,
        output_dir / "du_histogram.pdf",
    )
    assert output_path.exists()
    assert output_path.stat().st_size > 0

    bar_output_path = write_du_bar_chart(
        data_2026,
        output_dir / "du_bar_chart.pdf",
    )
    assert bar_output_path.exists()
    assert bar_output_path.stat().st_size > 0

    import matplotlib.pyplot as plt

    plt.close(figure)
    plt.close(figure_2026)
    plt.close(bar_figure)

print("Degree-of-unsaturation histogram verification completed")
