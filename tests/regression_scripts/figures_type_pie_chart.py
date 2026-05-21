from pathlib import Path
from tempfile import TemporaryDirectory

from baseline import load_production_baseline

from astromol.census import CensusView
from astromol.database import Database
from astromol.figures import (
    molecule_type_data,
    plot_type_pie_chart,
    write_type_pie_chart,
)


db = Database()
counts_2026 = load_production_baseline()["regression_counts"][
    "figures_molecule_type_2026"
]
view_2021 = CensusView.for_census(db, "2021")
data_2021 = molecule_type_data(view_2021)

assert data_2021.molecule_count == 240
assert data_2021.counts == {
    "neutral": 204,
    "radical": 54,
    "cation": 30,
    "cyclic": 19,
    "anion": 6,
    "fullerene": 3,
    "pah": 3,
}
assert {
    key: round(100 * fraction, 1)
    for key, fraction in data_2021.fractions.items()
} == {
    "neutral": 85.0,
    "radical": 22.5,
    "cation": 12.5,
    "cyclic": 7.9,
    "anion": 2.5,
    "fullerene": 1.2,
    "pah": 1.2,
}

view_2026 = CensusView.for_census(db, "2026")
data_2026 = molecule_type_data(view_2026)
assert data_2026.molecule_count == counts_2026["molecule_count"]
assert data_2026.counts == counts_2026["counts"]

with TemporaryDirectory() as tmp:
    output_dir = Path(tmp)
    figure, ax = plot_type_pie_chart(data_2021)
    displayed_text = [text.get_text() for text in ax.texts if text.get_text()]
    labels = displayed_text[:7]
    percents = displayed_text[7:]
    assert labels == [
        "Neutrals",
        "Radicals",
        "Cations",
        "Cyclics",
        "Anions",
        "Fullerenes",
        "PAHs",
    ]
    assert percents == [
        "85.0%",
        "22.5%",
        "12.5%",
        "7.9%",
        "2.5%",
        "1.2%",
        "1.2%",
    ]
    assert ax.axison is False
    plt_figure, plt_ax = plot_type_pie_chart(data_2026)
    displayed_2026 = [
        text.get_text()
        for text in plt_ax.texts
        if text.get_text()
    ]
    assert displayed_2026 == counts_2026["ring_text"]

    output_path = write_type_pie_chart(
        data_2021,
        output_dir / "type_pie_chart.pdf",
    )
    assert output_path.exists()
    assert output_path.stat().st_size > 0

    import matplotlib.pyplot as plt

    plt.close(figure)
    plt.close(plt_figure)

print("Molecule-type pie chart verification completed")
