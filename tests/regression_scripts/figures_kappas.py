from pathlib import Path
from tempfile import TemporaryDirectory

from baseline import load_production_baseline

from astromol.census import CensusView
from astromol.database import Database
from astromol.figures import (
    kappa_histogram_data,
    plot_kappa_histogram,
    write_kappa_histogram,
)


db = Database()
counts_current = load_production_baseline()["regression_counts"][
    "figures_kappas_current"
]
view_2021 = CensusView.for_census(db, "2021")
data_2021 = kappa_histogram_data(view_2021)

assert data_2021.molecule_count == 221
assert data_2021.min_kappa == -1.0
assert data_2021.max_kappa == 1.0
assert data_2021.histogram_counts().max() == 161
assert data_2021.histogram_counts()[:10].tolist() == [
    161,
    10,
    10,
    4,
    2,
    3,
    3,
    2,
    0,
    1,
]

view_current = CensusView.current(db, end_year=load_production_baseline()["analysis_end_year"])
data_current = kappa_histogram_data(view_current)
assert data_current.molecule_count == counts_current["molecule_count"]
assert data_current.min_kappa == -1.0
assert data_current.max_kappa == 1.0
assert data_current.histogram_counts().max() == counts_current["histogram_max"]

with TemporaryDirectory() as tmp:
    output_dir = Path(tmp)

    figure, ax = plot_kappa_histogram(data_2021)
    assert tuple(figure.get_size_inches()) == (8.0, 4.8)
    assert ax.get_title() == ""
    assert ax.get_xlabel() == "$\\kappa$"
    assert ax.get_ylabel() == "# Molecules"
    assert ax.get_yscale() == "log"
    assert list(ax.get_xticks()) == [-1.0, -0.5, 0.0, 0.5, 1.0]
    assert ax.xaxis.get_ticks_position() == "default"
    assert ax.yaxis.get_ticks_position() == "default"
    assert ax.spines["top"].get_visible()
    assert ax.spines["right"].get_visible()
    assert ax.spines["top"].get_edgecolor() == (0.0, 0.0, 0.0, 1.0)
    assert ax.spines["right"].get_edgecolor() == (0.0, 0.0, 0.0, 1.0)
    assert {text.get_text() for text in ax.texts} >= {
        "prolate",
        "asymmetric",
        "oblate",
    }

    output_path = write_kappa_histogram(
        data_2021,
        output_dir / "kappas.pdf",
    )
    assert output_path.exists()
    assert output_path.stat().st_size > 0

    import matplotlib.pyplot as plt

    plt.close(figure)

print("Kappa histogram verification completed")
