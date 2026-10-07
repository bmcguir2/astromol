from rendering_assertions import boxplot_labels
from pathlib import Path
from tempfile import TemporaryDirectory

from baseline import load_production_baseline

from astromol.census import CensusView
from astromol.database import Database
from astromol.figures import (
    mass_by_source_type_data,
    plot_mass_by_source_type,
    plot_mass_by_source_type_boxplot,
    write_mass_by_source_type,
    write_mass_by_source_type_boxplot,
)


db = Database()
counts_current = load_production_baseline()["regression_counts"][
    "figures_mass_by_source_type_current"
]
view_2021 = CensusView.for_census(db, "2021")
data_2021 = mass_by_source_type_data(view_2021)

assert data_2021.molecule_count == 217
assert data_2021.counts == {
    "carbon_star": 58,
    "dark_cloud": 67,
    "diffuse_cloud": 24,
    "sfr": 87,
}
assert tuple(round(value, 3) for value in data_2021.mass_range) == (2.016, 153.058)

view_current = CensusView.current(db, end_year=load_production_baseline()["analysis_end_year"])
data_current = mass_by_source_type_data(view_current)
assert data_current.molecule_count == counts_current["molecule_count"]
assert data_current.counts == counts_current["counts"]
assert [
    round(value, 3) for value in data_current.mass_range
] == counts_current["mass_range"]

with TemporaryDirectory() as tmp:
    output_dir = Path(tmp)

    figure, ax = plot_mass_by_source_type(
        data_2021,
        source_label_overrides={"diffuse_cloud": "LOS Cloud"},
    )
    displayed_text = [
        text.get_text()
        for text in ax.texts
        if text.get_text()
    ]
    assert displayed_text == [
        "Source Types",
        "SFR",
        "Carbon Star",
        "Dark Cloud",
        "LOS Cloud",
    ]
    assert ax.get_xlabel() == "Molecular Mass (amu)"
    assert ax.get_ylabel() == "Probability Density Estimate"

    output_path = write_mass_by_source_type(
        data_2021,
        output_dir / "mass_by_source_type_kde.pdf",
        source_label_overrides={"diffuse_cloud": "LOS Cloud"},
    )
    assert output_path.exists()
    assert output_path.stat().st_size > 0

    box_figure, box_ax = plot_mass_by_source_type_boxplot(data_current)
    assert box_ax.get_xlabel() == "Molecular Mass (amu)"
    assert [
        label.get_text()
        for label in box_ax.get_yticklabels()
    ] == [
        "Dark Cloud",
        "Carbon Star",
        "SFR",
        "Diffuse Cloud",
    ]
    assert [
        text.get_text()
        for text in box_ax.texts
        if text.get_text()
    ] == boxplot_labels(data_current, 'masses')

    box_output_path = write_mass_by_source_type_boxplot(
        data_current,
        output_dir / "mass_by_source_type_boxplot.pdf",
    )
    assert box_output_path.exists()
    assert box_output_path.stat().st_size > 0

    import matplotlib.pyplot as plt

    plt.close(figure)
    plt.close(box_figure)

print("Mass-by-source-type figure verification completed")
