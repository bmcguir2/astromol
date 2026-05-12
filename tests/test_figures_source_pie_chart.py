from pathlib import Path
from tempfile import TemporaryDirectory

from astromol.census import CensusView
from astromol.database import Database
from astromol.figures import (
    plot_source_pie_chart,
    source_type_data,
    write_source_pie_chart,
)


db = Database()
view_2021 = CensusView.for_census(db, "2021")
data_2021 = source_type_data(view_2021)

assert data_2021.molecule_count == 240
assert data_2021.counts == {
    "sfr": 87,
    "dark_cloud": 67,
    "carbon_star": 58,
    "other": 32,
    "diffuse_cloud": 24,
}
assert {
    key: round(100 * fraction, 1)
    for key, fraction in data_2021.fractions.items()
} == {
    "sfr": 36.2,
    "dark_cloud": 27.9,
    "carbon_star": 24.2,
    "other": 13.3,
    "diffuse_cloud": 10.0,
}

view_2026 = CensusView.for_census(db, "2026")
data_2026 = source_type_data(view_2026)
assert data_2026.molecule_count == 325
assert data_2026.counts == {
    "sfr": 92,
    "dark_cloud": 117,
    "carbon_star": 71,
    "other": 50,
    "diffuse_cloud": 24,
}

with TemporaryDirectory() as tmp:
    output_dir = Path(tmp)
    figure, ax = plot_source_pie_chart(
        data_2021,
        diffuse_cloud_label="LOS Cloud",
    )
    displayed_text = [text.get_text() for text in ax.texts if text.get_text()]
    assert displayed_text == [
        "SFR",
        "36.2%",
        "Dark Cloud",
        "27.9%",
        "Carbon Star",
        "24.2%",
        "Other",
        "13.3%",
        "LOS Cloud",
        "10.0%",
    ]
    assert ax.axison is False
    production_figure, production_ax = plot_source_pie_chart(data_2026)
    displayed_2026 = [
        text.get_text()
        for text in production_ax.texts
        if text.get_text()
    ]
    assert displayed_2026 == [
        "Dark Cloud",
        "36.0%",
        "SFR",
        "28.3%",
        "Carbon Star",
        "21.8%",
        "Other",
        "15.4%",
        "Diffuse Cloud",
        "7.4%",
    ]

    output_path = write_source_pie_chart(
        data_2021,
        output_dir / "source_pie_chart.pdf",
        diffuse_cloud_label="LOS Cloud",
    )
    assert output_path.exists()
    assert output_path.stat().st_size > 0

    import matplotlib.pyplot as plt

    plt.close(figure)
    plt.close(production_figure)

print("Source-type pie chart verification completed")
