from pathlib import Path
from tempfile import TemporaryDirectory

from baseline import load_production_baseline

from astromol.census import CensusView
from astromol.database import Database
from astromol.figures import (
    MIT_RED,
    plot_scopes_by_year,
    scopes_by_year_data,
    write_scopes_by_year_plot,
)


def rows_by_label(data):
    """Return scope rows keyed by plotted facility label."""
    return {series.label: series for series in data.series}


db = Database()
counts_2026 = load_production_baseline()["regression_counts"][
    "figures_scopes_by_year_2026"
]
view_2021 = CensusView.for_census(db, "2021")
view_2026 = CensusView.for_census(db, "2026")

data_2021 = scopes_by_year_data(view_2021)
rows_2021 = rows_by_label(data_2021)

assert data_2021.end_year == 2021
assert len(data_2021.series) == 7
assert rows_2021["NRAO 36-ft"].final_count == 33
assert rows_2021["IRAM 30-m"].final_count == 64
assert rows_2021["GBT 100-m"].final_count == 28
assert rows_2021["Nobeyama 45-m"].final_count == 15
assert rows_2021["NRAO/ARO 12-m"].final_count == 27
assert rows_2021["Yebes 40-m"].final_count == 19
assert rows_2021["NRAO 140-ft"].final_count == 13

assert round(rows_2021["NRAO 36-ft"].rate, 1) == 2.2
assert round(rows_2021["IRAM 30-m"].rate, 1) == 1.5
assert round(rows_2021["GBT 100-m"].rate, 1) == 1.1
assert round(rows_2021["Nobeyama 45-m"].rate, 1) == 1.1
assert round(rows_2021["NRAO/ARO 12-m"].rate, 1) == 0.8
assert round(rows_2021["Yebes 40-m"].rate, 1) == 0.7
assert round(rows_2021["NRAO 140-ft"].rate, 1) == 0.5

assert rows_2021["NRAO 36-ft"].fit_stop_year == 1985
assert rows_2021["Nobeyama 45-m"].fit_stop_year == 1997
assert rows_2021["NRAO 140-ft"].fit_stop_year == 1993

data_2026 = scopes_by_year_data(view_2026)
rows_2026 = rows_by_label(data_2026)
assert data_2026.end_year == 2026
assert rows_2026["Yebes 40-m"].final_count == counts_2026["yebes_final_count"]
assert round(rows_2026["Yebes 40-m"].rate, 1) == counts_2026["yebes_rate_rounded"]
assert "ALMA" in rows_2026

with TemporaryDirectory() as tmp:
    output_dir = Path(tmp)
    figure, ax = plot_scopes_by_year(data_2021)
    assert ax.get_title() == ""
    assert ax.get_xlabel() == "Year"
    assert ax.get_ylabel() == "Cumulative Number of Detected Molecules"
    assert ax.get_xlim() == (1965.0, 2023.0)
    assert all(spine.get_visible() for spine in ax.spines.values())
    assert all(
        spine.get_edgecolor() == (0.0, 0.0, 0.0, 1.0)
        for spine in ax.spines.values()
    )
    assert any(text.get_text() == "NRAO 36-ft" for text in ax.texts)
    assert any(text.get_text() == "2.2/yr" for text in ax.texts)
    assert any(text.get_text() == "(1967 - 1985)" for text in ax.texts)

    modern_figure, modern_ax = plot_scopes_by_year(data_2026, style="modern")
    yebes_lines = [
        line
        for line in modern_ax.lines
        if line.get_color().lower() == MIT_RED.lower()
    ]
    assert len(yebes_lines) == 1
    assert yebes_lines[0].get_linewidth() == 4.0
    assert yebes_lines[0].get_xdata()[0] >= 2007

    dormant_lines = [
        line
        for line in modern_ax.lines
        if line.get_linestyle() != "-"
    ]
    assert dormant_lines
    assert all(line.get_alpha() == 0.50 for line in dormant_lines)

    output_path = write_scopes_by_year_plot(
        data_2021,
        output_dir / "scopes_by_year.pdf",
    )
    assert output_path.exists()
    assert output_path.stat().st_size > 0

    modern_output_path = write_scopes_by_year_plot(
        data_2026,
        output_dir / "scopes_by_year_modern.pdf",
        style="modern",
    )
    assert modern_output_path.exists()
    assert modern_output_path.stat().st_size > 0

    import matplotlib.pyplot as plt

    plt.close(figure)
    plt.close(modern_figure)

print("Scopes-by-year figure verification passed")
