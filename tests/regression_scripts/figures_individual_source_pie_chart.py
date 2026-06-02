from pathlib import Path
from tempfile import TemporaryDirectory

from baseline import load_production_baseline

from astromol.census import CensusView
from astromol.database import Database
from astromol.figures import (
    individual_source_data,
    plot_individual_source_pie_chart,
    write_individual_source_pie_chart,
)


def ring_text(data):
    categories = tuple(
        category
        for _, category in sorted(
            enumerate(data.categories),
            key=lambda indexed_category: (
                -indexed_category[1].count,
                indexed_category[0],
            ),
        )
    )
    text = []
    for category in categories:
        text.extend([category.label, f"{category.percent:.1f}%"])
    return text


db = Database()
counts_2026 = load_production_baseline()["regression_counts"][
    "figures_individual_source_2026"
]
view_2021 = CensusView.for_census(db, "2021")
data_2021 = individual_source_data(view_2021)

assert set(data_2021.counts) == {
    "other",
    "sgr_b2",
    "tmc1",
    "irc10216",
    "orion",
    "g0693",
}
assert next(
    category for category in data_2021.categories if category.key == "g0693"
).color == "black"

view_2026 = CensusView.for_census(db, "2026")
data_2026 = individual_source_data(view_2026)
assert data_2026.molecule_count == counts_2026["molecule_count"]
assert data_2026.counts == counts_2026["counts"]
assert next(
    category for category in data_2026.categories if category.key == "g0693"
).color == "black"

with TemporaryDirectory() as tmp:
    output_dir = Path(tmp)
    figure, ax = plot_individual_source_pie_chart(data_2021)
    displayed_text = [text.get_text() for text in ax.texts if text.get_text()]
    assert displayed_text == ring_text(data_2021)
    g0693_label = next(text for text in ax.texts if text.get_text() == "G+0.693")
    assert g0693_label.get_color() == "black"
    assert ax.axison is False

    production_figure, production_ax = plot_individual_source_pie_chart(data_2026)
    displayed_2026 = [
        text.get_text()
        for text in production_ax.texts
        if text.get_text()
    ]
    assert displayed_2026 == counts_2026["ring_text"]
    production_g0693_label = next(
        text for text in production_ax.texts if text.get_text() == "G+0.693"
    )
    assert production_g0693_label.get_color() == "black"

    output_path = write_individual_source_pie_chart(
        data_2021,
        output_dir / "indiv_source_pie_chart.pdf",
    )
    assert output_path.exists()
    assert output_path.stat().st_size > 0

    import matplotlib.pyplot as plt

    plt.close(figure)
    plt.close(production_figure)

print("Individual-source pie chart verification completed")
