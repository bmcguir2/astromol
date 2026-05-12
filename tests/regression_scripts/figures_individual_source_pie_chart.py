from pathlib import Path
from tempfile import TemporaryDirectory

from astromol.census import CensusView
from astromol.database import Database
from astromol.figures import (
    individual_source_data,
    plot_individual_source_pie_chart,
    write_individual_source_pie_chart,
)


db = Database()
view_2021 = CensusView.for_census(db, "2021")
data_2021 = individual_source_data(view_2021)

assert data_2021.molecule_count == 240
assert data_2021.counts == {
    "other": 128,
    "sgr_b2": 68,
    "tmc1": 57,
    "irc10216": 55,
    "orion": 24,
}
assert {
    key: round(100 * fraction, 1)
    for key, fraction in data_2021.fractions.items()
} == {
    "other": 53.3,
    "sgr_b2": 28.3,
    "tmc1": 23.8,
    "irc10216": 22.9,
    "orion": 10.0,
}

view_2026 = CensusView.for_census(db, "2026")
data_2026 = individual_source_data(view_2026)
assert data_2026.molecule_count == 325
assert data_2026.counts == {
    "other": 155,
    "sgr_b2": 70,
    "tmc1": 106,
    "irc10216": 68,
    "orion": 25,
}

with TemporaryDirectory() as tmp:
    output_dir = Path(tmp)
    figure, ax = plot_individual_source_pie_chart(data_2021)
    displayed_text = [text.get_text() for text in ax.texts if text.get_text()]
    assert displayed_text == [
        "Other",
        "53.3%",
        "Sgr B2",
        "28.3%",
        "TMC-1",
        "23.8%",
        "IRC+10216",
        "22.9%",
        "Orion",
        "10.0%",
    ]
    assert ax.axison is False

    production_figure, production_ax = plot_individual_source_pie_chart(data_2026)
    displayed_2026 = [
        text.get_text()
        for text in production_ax.texts
        if text.get_text()
    ]
    assert displayed_2026 == [
        "Other",
        "47.7%",
        "TMC-1",
        "32.6%",
        "Sgr B2",
        "21.5%",
        "IRC+10216",
        "20.9%",
        "Orion",
        "7.7%",
    ]

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
