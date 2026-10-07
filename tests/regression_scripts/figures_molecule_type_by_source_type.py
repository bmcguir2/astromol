from pathlib import Path
from tempfile import TemporaryDirectory
from collections import Counter

from baseline import load_production_baseline

from astromol.census import CensusView
from astromol.database import Database
from astromol.figures import (
    molecule_type_by_source_type_data,
    plot_molecule_type_by_source_enrichment_matrix,
    plot_molecule_type_by_source_type,
    write_molecule_type_by_source_enrichment_matrix,
    write_molecule_type_by_source_type,
)


db = Database()
counts_current = load_production_baseline()["regression_counts"][
    "figures_molecule_type_by_source_type_current"
]
view_2021 = CensusView.for_census(db, "2021")
data_2021 = molecule_type_by_source_type_data(view_2021)

assert data_2021.molecule_count == 240
assert data_2021.counts == {
    "carbon_star": {
        "anion": 6,
        "cation": 0,
        "cyclic": 3,
        "neutral": 52,
        "radical": 19,
    },
    "dark_cloud": {
        "anion": 2,
        "cation": 9,
        "cyclic": 10,
        "neutral": 56,
        "radical": 17,
    },
    "diffuse_cloud": {
        "anion": 0,
        "cation": 6,
        "cyclic": 0,
        "neutral": 18,
        "radical": 9,
    },
    "sfr": {
        "anion": 0,
        "cation": 13,
        "cyclic": 4,
        "neutral": 74,
        "radical": 7,
    },
}
assert data_2021.source_counts == {
    "carbon_star": 58,
    "dark_cloud": 67,
    "diffuse_cloud": 24,
    "sfr": 87,
}
assert data_2021.overall_type_counts == {
    "anion": 6,
    "cation": 30,
    "cyclic": 19,
    "neutral": 204,
    "radical": 54,
}

view_current = CensusView.current(db, end_year=load_production_baseline()["analysis_end_year"])
data_current = molecule_type_by_source_type_data(view_current)
assert data_current.molecule_count == counts_current["molecule_count"]
assert data_current.counts == counts_current["counts"]
assert data_current.source_counts == counts_current["source_counts"]
assert data_current.overall_type_counts == counts_current["overall_type_counts"]

with TemporaryDirectory() as tmp:
    output_dir = Path(tmp)
    figure, axes = plot_molecule_type_by_source_type(
        data_2021,
        source_label_overrides={"diffuse_cloud": "LOS Clouds"},
    )
    axis_texts = [
        [text.get_text() for text in ax.texts if text.get_text()]
        for ax in axes.flat
    ]
    assert axis_texts == [
        ["6", "3", "52", "19", "Carbon Stars"],
        ["2", "9", "10", "56", "17", "Dark Clouds"],
        ["6", "18", "9", "LOS Clouds"],
        ["13", "4", "74", "7", "SFRs"],
    ]
    assert len(axes[0, 1].get_legend().texts) == 5

    output_path = write_molecule_type_by_source_type(
        data_2021,
        output_dir / "mol_type_by_source_type.pdf",
        source_label_overrides={"diffuse_cloud": "LOS Clouds"},
    )
    assert output_path.exists()
    assert output_path.stat().st_size > 0

    matrix_figure, matrix_ax = plot_molecule_type_by_source_enrichment_matrix(
        data_current,
    )
    matrix_text = [
        text.get_text()
        for text in matrix_ax.texts
        if text.get_text()
    ]
    assert matrix_text and all("nan" not in t.lower() and "inf" not in t.lower() for t in matrix_text)
    assert Counter(t for t in matrix_text if t.startswith("n=")) == Counter(
        f"n={count}" for counts in data_current.counts.values() for count in counts.values()
    )
    assert [
        label.get_text()
        for label in matrix_ax.get_xticklabels()
    ] == [
        "Neutral",
        "Radical",
        "Cation",
        "Anion",
        "Cyclic",
    ]

    matrix_output_path = write_molecule_type_by_source_enrichment_matrix(
        data_current,
        output_dir / "mol_type_by_source_enrichment_matrix.pdf",
    )
    assert matrix_output_path.exists()
    assert matrix_output_path.stat().st_size > 0

    import matplotlib.pyplot as plt

    plt.close(figure)
    plt.close(matrix_figure)

print("Molecule-type-by-source figure verification completed")
