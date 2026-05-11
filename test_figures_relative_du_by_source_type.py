from pathlib import Path
from tempfile import TemporaryDirectory

from astromol.census import CensusView
from astromol.database import Database
from astromol.figures import (
    plot_relative_du_by_source_type,
    plot_relative_du_by_source_type_boxplot,
    relative_du_by_source_type_data,
    write_relative_du_by_source_type,
    write_relative_du_by_source_type_boxplot,
)


db = Database()
view_2021 = CensusView.for_census(db, "2021")
data_2021 = relative_du_by_source_type_data(view_2021)

assert data_2021.molecule_count == 184
assert data_2021.counts == {
    "carbon_star": 27,
    "dark_cloud": 66,
    "diffuse_cloud": 24,
    "sfr": 85,
}
assert {
    key: (round(min(values), 3), round(max(values), 3))
    for key, values in data_2021.values_by_source.items()
} == {
    "carbon_star": (0.333, 1.0),
    "dark_cloud": (-0.333, 1.0),
    "diffuse_cloud": (-0.5, 1.0),
    "sfr": (-0.5, 1.0),
}

view_2026 = CensusView.for_census(db, "2026")
data_2026 = relative_du_by_source_type_data(view_2026)
assert data_2026.molecule_count == 239
assert data_2026.counts == {
    "carbon_star": 28,
    "dark_cloud": 116,
    "diffuse_cloud": 24,
    "sfr": 90,
}

with TemporaryDirectory() as tmp:
    output_dir = Path(tmp)
    figure, axes = plot_relative_du_by_source_type(
        data_2021,
        source_label_overrides={"diffuse_cloud": "LOS Cloud"},
    )
    axis_texts = [
        [text.get_text() for text in ax.texts if text.get_text()]
        for ax in axes
    ]
    assert axis_texts == [
        ["Carbon Star"],
        ["LOS Cloud"],
        ["Dark Cloud"],
        ["SFR"],
    ]
    assert axes[2].get_xlabel() == "Relative Degree of Unsaturation"
    assert axes[2].get_ylabel() == "Probability Density Estimate"

    output_path = write_relative_du_by_source_type(
        data_2021,
        output_dir / "relative_du_by_source_type_kde.pdf",
        source_label_overrides={"diffuse_cloud": "LOS Cloud"},
    )
    assert output_path.exists()
    assert output_path.stat().st_size > 0

    box_figure, box_ax = plot_relative_du_by_source_type_boxplot(data_2026)
    assert box_ax.get_xlabel() == "Relative Degree of Unsaturation"
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
    ] == [
        "n=115",
        "n=28",
        "n=87",
        "n=22",
    ]

    box_output_path = write_relative_du_by_source_type_boxplot(
        data_2026,
        output_dir / "relative_du_by_source_type_boxplot.pdf",
    )
    assert box_output_path.exists()
    assert box_output_path.stat().st_size > 0

    import matplotlib.pyplot as plt

    plt.close(figure)
    plt.close(box_figure)

print("Relative-DU-by-source-type figure verification completed")
