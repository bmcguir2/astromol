import os
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from baseline import load_production_baseline

from astromol.census import CensusView
from astromol.database import Database
from astromol.figures import (
    cumulative_detection_data,
    plot_cumulative_detections,
    write_cumulative_detections_plot,
)


def count_at_year(data, year):
    """Return the cumulative count at a specific year."""
    index = int(np.argwhere(data.years == year)[0][0])
    return int(data.counts[index])


db = Database()
counts_2026 = load_production_baseline()["regression_counts"][
    "figures_cumulative_detections_2026"
]
view_2021 = CensusView.for_census(db, "2021")
view_2026 = CensusView.for_census(db, "2026")
view_current = CensusView.current(db)

data_2021 = cumulative_detection_data(view_2021)
assert data_2021.start_year == 1937
assert data_2021.end_year == 2021
assert data_2021.total == 240
assert len(data_2021.first_detection_years) == 240
assert len(data_2021.years) == 85
assert data_2021.show_facility_markers is True

assert count_at_year(data_2021, 1937) == 1
assert count_at_year(data_2021, 1940) == 2
assert count_at_year(data_2021, 1941) == 3
assert count_at_year(data_2021, 1963) == 4
assert count_at_year(data_2021, 1968) == 5
assert count_at_year(data_2021, 2004) == 130
assert count_at_year(data_2021, 2005) == 130
assert count_at_year(data_2021, 2019) == 213
assert count_at_year(data_2021, 2020) == 219
assert count_at_year(data_2021, 2021) == 240

trends_2021 = {trend.label: trend.slope for trend in data_2021.trends}
assert set(trends_2021) == {"1968-2005", "2005-2021"}
assert np.isclose(trends_2021["1968-2005"], 3.4569427727322477)
assert np.isclose(trends_2021["2005-2021"], 5.9779411764705905)

data_2026 = cumulative_detection_data(view_2026)
assert data_2026.start_year == 1937
assert data_2026.end_year == 2026
assert data_2026.total == counts_2026["total"]
assert len(data_2026.first_detection_years) == counts_2026["first_detection_years"]
assert data_2026.show_facility_markers is False
assert count_at_year(data_2026, 2019) == counts_2026["counts_by_year"]["2019"]
assert count_at_year(data_2026, 2020) == counts_2026["counts_by_year"]["2020"]
assert count_at_year(data_2026, 2021) == counts_2026["counts_by_year"]["2021"]
assert count_at_year(data_2026, 2024) == counts_2026["counts_by_year"]["2024"]
assert count_at_year(data_2026, 2026) == counts_2026["counts_by_year"]["2026"]

trends_2026 = {trend.label: trend.slope for trend in data_2026.trends}
assert set(trends_2026) == {"1968-2005", "2005-2021", "2021-2026"}
assert np.isclose(trends_2026["1968-2005"], 3.4569427727322477)
assert np.isclose(trends_2026["2005-2021"], 6.394607843137259)
assert np.isclose(
    trends_2026["2021-2026"],
    counts_2026["trend_slopes"]["2021-2026"],
)

data_current = cumulative_detection_data(view_current)
trends_current = {trend.label: trend.slope for trend in data_current.trends}
assert data_current.show_facility_markers is False
assert set(trends_current) == {"1968-2005", "2005-2021", "2021-Present"}
assert np.isclose(
    trends_current["2021-Present"],
    counts_2026["trend_slopes"]["2021-2026"],
)

with TemporaryDirectory() as tmp:
    output_dir = Path(tmp)
    os.environ.setdefault("MPLCONFIGDIR", str(output_dir / "mplconfig"))
    os.environ.setdefault("XDG_CACHE_HOME", str(output_dir / "cache"))

    figure, axes = plot_cumulative_detections(data_2021)
    assert axes.get_ylim() == (0.0, 270.0)
    assert axes.texts[0].get_position() == (0.05, 0.95)
    assert axes.texts[1].get_position() == (0.95, 0.95)

    output_path = write_cumulative_detections_plot(
        data_2021,
        output_dir / "cumulative_detections.pdf",
    )
    assert output_path.exists()
    assert output_path.stat().st_size > 0

    import matplotlib.pyplot as plt

    plt.close(figure)

print("Cumulative detections figure data verification passed")
