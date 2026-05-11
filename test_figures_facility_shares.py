from pathlib import Path
from tempfile import TemporaryDirectory

from astromol.census import CensusView
from astromol.database import Database
from astromol.figures import (
    FACILITY_SHARE_INACTIVE,
    LEGACY_2021_FACILITY_SHARE_NICKS,
    facility_share_data,
    plot_facility_share_bars,
    plot_facility_shares,
    write_facility_share_bars_plot,
    write_facility_shares_plot,
)


def rows_by_nick(data):
    """Return facility-share rows keyed by telescope nick."""
    return {facility.nick: facility for facility in data.facilities}


db = Database()
view_2021 = CensusView.for_census(db, "2021")
view_2026 = CensusView.for_census(db, "2026")

data_2021 = facility_share_data(view_2021)
assert data_2021.end_year == 2021
assert tuple(facility.nick for facility in data_2021.facilities) == (
    "NRAO36",
    "IRAM30",
    "GBT",
    "Herschel",
    "Yebes40",
    "NRAOARO12",
    "Bell7m",
    "ALMA",
    "NRAO140",
)

rows_2021 = rows_by_nick(data_2021)
assert rows_2021["NRAO36"].detection_count == 33
assert rows_2021["NRAO36"].total_window_detections == 56
assert rows_2021["NRAO36"].percent == 58
assert rows_2021["IRAM30"].percent == 34
assert rows_2021["GBT"].percent == 24
assert rows_2021["Herschel"].percent == 21
assert rows_2021["Yebes40"].percent == 18
assert rows_2021["NRAOARO12"].percent == 14
assert rows_2021["Bell7m"].percent == 13
assert rows_2021["ALMA"].percent == 9
assert rows_2021["NRAO140"].percent == 8
assert rows_2021["Herschel"].active_at_view is False
assert rows_2021["ALMA"].active_at_view is True

modern_selection_2021 = facility_share_data(
    view_2021,
    use_legacy_2021_selection=False,
)
modern_selection_nicks = {facility.nick for facility in modern_selection_2021.facilities}
assert "Nobeyama45" in modern_selection_nicks
assert "Herschel" not in modern_selection_nicks
assert set(LEGACY_2021_FACILITY_SHARE_NICKS) - modern_selection_nicks == {"Herschel"}

data_2026 = facility_share_data(view_2026)
assert data_2026.end_year == 2026
assert "Nobeyama45" in {facility.nick for facility in data_2026.facilities}

with TemporaryDirectory() as tmp:
    output_dir = Path(tmp)
    figure, axes = plot_facility_shares(data_2021)
    assert axes.shape == (3, 3)
    assert len(figure.axes) == 9
    assert axes[0, 0].texts[0].get_text() == "58%"
    assert axes[0, 0].texts[0].get_fontweight() == "bold"
    assert axes[0, 0].texts[0].get_fontsize() == 15
    assert axes[0, 0].texts[2].get_text() == "NRAO 36-ft"
    assert axes[0, 0].texts[2].get_fontsize() == 17
    assert axes[0, 0].texts[2].get_fontweight() == "normal"
    assert axes[0, 0].texts[2].get_position() == (0.5, 1.02)
    assert axes[0, 0].texts[3].get_text() == "1967 - 1984"
    assert axes[0, 0].texts[3].get_fontsize() == 14
    assert axes[0, 0].texts[3].get_fontweight() == "normal"
    assert axes[0, 0].texts[3].get_position() == (0.5, 0.925)
    assert axes[0, 0].patches[0].r == 0.9
    assert FACILITY_SHARE_INACTIVE == "#F87070"

    output_path = write_facility_shares_plot(
        data_2021,
        output_dir / "facility_shares.pdf",
    )
    assert output_path.exists()
    assert output_path.stat().st_size > 0

    bar_figure, bar_ax = plot_facility_share_bars(data_2021)
    assert bar_ax.get_title() == ""
    assert bar_ax.get_xlabel() == "Contribution share during facility lifetime (%)"
    assert all(not line.get_visible() for line in bar_ax.xaxis.get_gridlines())
    bar_texts = {text.get_text(): text for text in bar_ax.texts}
    assert "NRAO 36-ft" in bar_texts
    assert bar_texts["NRAO 36-ft"].get_fontweight() == "bold"
    assert "1967 - 1984" in bar_texts
    assert bar_texts["1967 - 1984"].get_fontstyle() == "italic"
    assert "58% (33/56)" in bar_texts
    assert bar_texts["58% (33/56)"].get_fontweight() == "bold"
    assert all(spine.get_visible() for spine in bar_ax.spines.values())
    assert all(
        spine.get_edgecolor() == (0.0, 0.0, 0.0, 1.0)
        for spine in bar_ax.spines.values()
    )
    legend = bar_ax.get_legend()
    assert legend is not None
    assert [text.get_text() for text in legend.get_texts()] == [
        "Active",
        "Decommissioned",
    ]
    assert legend.handleheight == 2.4

    bar_output_path = write_facility_share_bars_plot(
        data_2021,
        output_dir / "facility_share_bars.pdf",
    )
    assert bar_output_path.exists()
    assert bar_output_path.stat().st_size > 0

    import matplotlib.pyplot as plt

    plt.close(figure)
    plt.close(bar_figure)

print("Facility-shares figure verification passed")
