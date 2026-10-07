from rendering_assertions import boxplot_labels
from pathlib import Path
from tempfile import TemporaryDirectory

from baseline import load_production_baseline

from astromol.census import CensusView
from astromol.database import Database
from astromol.figures import (
    du_by_source_type_data,
    plot_du_by_source_type,
    plot_du_by_source_type_boxplot,
    write_du_by_source_type_boxplot,
    write_du_by_source_type,
)


db = Database()
counts_current = load_production_baseline()["regression_counts"][
    "figures_du_by_source_type_current"
]
view_2021 = CensusView.for_census(db, "2021")
data_2021 = du_by_source_type_data(view_2021)

assert data_2021.molecule_count == 184
assert data_2021.counts == {
    "carbon_star": 27,
    "dark_cloud": 66,
    "diffuse_cloud": 24,
    "sfr": 85,
}

view_current = CensusView.current(db, end_year=load_production_baseline()["analysis_end_year"])
data_current = du_by_source_type_data(view_current)
assert data_current.molecule_count == counts_current["molecule_count"]
assert data_current.counts == counts_current["counts"]

with TemporaryDirectory() as tmp:
    output_dir = Path(tmp)
    figure, ax = plot_du_by_source_type(
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
        "85",
        "27",
        "66",
        "24",
    ]
    assert ax.get_xlabel() == "Degree of Unsaturation"
    assert ax.get_ylabel() == "Probability Density Estimate"

    output_path = write_du_by_source_type(
        data_2021,
        output_dir / "du_by_source_type_kde.pdf",
        source_label_overrides={"diffuse_cloud": "LOS Cloud"},
    )
    assert output_path.exists()
    assert output_path.stat().st_size > 0

    box_figure, box_ax = plot_du_by_source_type_boxplot(data_current)
    assert box_ax.get_xlabel() == "Degree of Unsaturation"
    assert box_ax.get_xlim()[0] < 0
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
    ] == boxplot_labels(data_current, 'values')

    box_output_path = write_du_by_source_type_boxplot(
        data_current,
        output_dir / "du_by_source_type_boxplot.pdf",
    )
    assert box_output_path.exists()
    assert box_output_path.stat().st_size > 0

    import matplotlib.pyplot as plt

    plt.close(figure)
    plt.close(box_figure)

print("DU-by-source-type figure verification completed")
