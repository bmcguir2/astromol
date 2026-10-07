from rendering_assertions import ring_text as expected_ring_text
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


LABEL_Y_POSITIONS = (0.110, 0.160, 0.205, 0.255, 0.305, 0.355)
PERCENT_Y_POSITIONS = (0.870, 0.825, 0.775, 0.725, 0.680, 0.635)


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


def assert_ring_text_positions(ax, expected_text):
    label_text = expected_text[::2]
    percent_text = expected_text[1::2]
    labels = [text for text in ax.texts if text.get_text() in label_text]
    percentages = [text for text in ax.texts if text.get_text() in percent_text]

    assert [text.get_text() for text in labels] == label_text
    assert [round(text.get_position()[1], 3) for text in labels] == list(
        LABEL_Y_POSITIONS
    )
    assert [text.get_text() for text in percentages] == percent_text
    assert [round(text.get_position()[1], 3) for text in percentages] == list(
        PERCENT_Y_POSITIONS
    )


db = Database()
counts_current = load_production_baseline()["regression_counts"][
    "figures_individual_source_current"
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

view_current = CensusView.current(db, end_year=load_production_baseline()["analysis_end_year"])
data_current = individual_source_data(view_current)
assert data_current.molecule_count == counts_current["molecule_count"]
assert data_current.counts == counts_current["counts"]
assert next(
    category for category in data_current.categories if category.key == "g0693"
).color == "black"

with TemporaryDirectory() as tmp:
    output_dir = Path(tmp)
    figure, ax = plot_individual_source_pie_chart(data_2021)
    displayed_text = [text.get_text() for text in ax.texts if text.get_text()]
    expected_2021_text = ring_text(data_2021)
    assert displayed_text == expected_2021_text
    assert_ring_text_positions(ax, expected_2021_text)
    g0693_label = next(text for text in ax.texts if text.get_text() == "G+0.693")
    assert g0693_label.get_color() == "black"
    assert ax.axison is False

    production_figure, production_ax = plot_individual_source_pie_chart(data_current)
    displayed_current = [
        text.get_text()
        for text in production_ax.texts
        if text.get_text()
    ]
    assert displayed_current == expected_ring_text(data_current, types=False)
    assert_ring_text_positions(production_ax, expected_ring_text(data_current, types=False))
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
