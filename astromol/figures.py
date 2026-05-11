"""Figure data builders and plotting helpers for census outputs."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from datetime import date
from pathlib import Path

import numpy as np
from molmass import ELEMENTS

from .census import CensusView
from .models import Detection


FIGURE_SIZE = (10, 8)
FIGURE_AXES_BOUNDS = (0.14, 0.13, 0.82, 0.82)
FIGURE_TEXT_SIZE = 24
FIGURE_LEGEND_TEXT_SIZE = 16
FIGURE_TICK_LENGTH = 15
FIGURE_TICK_WIDTH = 1

ASTROMOL_BLUE = "dodgerblue"
MIT_RED = "#750014"
NRAO_BLUE = "#0A1589"
MIT_GRAY = "#8B959E"

DU_HISTOGRAM_ALLOWED_ELEMENTS = frozenset({"H", "D", "N", "C", "Cl", "F", "S", "O"})
DU_HISTOGRAM_BINS = tuple(np.arange(-0.25, 12.5, 0.5))
TYPE_PIE_SPECS = (
    {
        "key": "neutral",
        "label": "Neutral",
        "plural_label": "Neutrals",
        "attribute": "neutral",
        "color": ASTROMOL_BLUE,
        "percent_color": "white",
    },
    {
        "key": "radical",
        "label": "Radical",
        "plural_label": "Radicals",
        "attribute": "radical",
        "color": "darkorange",
        "percent_color": "darkorange",
    },
    {
        "key": "cation",
        "label": "Cation",
        "plural_label": "Cations",
        "attribute": "cation",
        "color": "forestgreen",
        "percent_color": "forestgreen",
    },
    {
        "key": "cyclic",
        "label": "Cyclic",
        "plural_label": "Cyclics",
        "attribute": "cyclic",
        "color": "violet",
        "percent_color": "violet",
    },
    {
        "key": "anion",
        "label": "Anion",
        "plural_label": "Anions",
        "attribute": "anion",
        "color": "red",
        "percent_color": "red",
    },
    {
        "key": "fullerene",
        "label": "Fullerene",
        "plural_label": "Fullerenes",
        "attribute": "fullerene",
        "color": "goldenrod",
        "percent_color": "goldenrod",
    },
    {
        "key": "pah",
        "label": "PAH",
        "plural_label": "PAHs",
        "attribute": "pah",
        "color": "royalblue",
        "percent_color": "royalblue",
    },
)
SOURCE_PIE_SPECS = (
    {
        "key": "sfr",
        "label": "SFR",
        "source_types": ("SFR",),
        "color": ASTROMOL_BLUE,
        "label_y": 0.110,
        "percent_y": 0.870,
    },
    {
        "key": "dark_cloud",
        "label": "Dark Cloud",
        "source_types": ("Dark Cloud",),
        "color": "forestgreen",
        "label_y": 0.160,
        "percent_y": 0.825,
    },
    {
        "key": "carbon_star",
        "label": "Carbon Star",
        "source_types": ("Carbon Star",),
        "color": "darkorange",
        "label_y": 0.205,
        "percent_y": 0.775,
    },
    {
        "key": "other",
        "label": "Other",
        "source_types": (),
        "color": "violet",
        "label_y": 0.255,
        "percent_y": 0.725,
    },
    {
        "key": "diffuse_cloud",
        "label": "Diffuse Cloud",
        "source_types": ("Diffuse Cloud",),
        "color": "red",
        "label_y": 0.305,
        "percent_y": 0.680,
    },
)
SOURCE_PIE_DIRECT_SOURCE_TYPES = {
    source_type
    for spec in SOURCE_PIE_SPECS
    for source_type in spec["source_types"]
}
INDIVIDUAL_SOURCE_PIE_SPECS = (
    {
        "key": "other",
        "label": "Other",
        "source_nicks": (),
        "color": ASTROMOL_BLUE,
    },
    {
        "key": "sgr_b2",
        "label": "Sgr B2",
        "source_nicks": ("SgrB2",),
        "color": "darkorange",
    },
    {
        "key": "tmc1",
        "label": "TMC-1",
        "source_nicks": ("TMC1",),
        "color": "violet",
    },
    {
        "key": "irc10216",
        "label": "IRC+10216",
        "source_nicks": ("IRC10216",),
        "color": "forestgreen",
    },
    {
        "key": "orion",
        "label": "Orion",
        "source_nicks": ("OrionKL",),
        "color": "red",
    },
)
INDIVIDUAL_SOURCE_NICK_TO_KEY = {
    source_nick: spec["key"]
    for spec in INDIVIDUAL_SOURCE_PIE_SPECS
    for source_nick in spec["source_nicks"]
}
MOLECULE_TYPE_BY_SOURCE_SOURCE_SPECS = (
    {
        "key": "carbon_star",
        "label": "Carbon Stars",
    },
    {
        "key": "dark_cloud",
        "label": "Dark Clouds",
    },
    {
        "key": "diffuse_cloud",
        "label": "Diffuse Clouds",
    },
    {
        "key": "sfr",
        "label": "SFRs",
    },
)
MOLECULE_TYPE_BY_SOURCE_TYPE_SPECS = (
    {
        "key": "anion",
        "label": "Anion",
        "attribute": "anion",
        "color": "darkorange",
    },
    {
        "key": "cation",
        "label": "Cation",
        "attribute": "cation",
        "color": "forestgreen",
    },
    {
        "key": "cyclic",
        "label": "Cyclic",
        "attribute": "cyclic",
        "color": "violet",
    },
    {
        "key": "neutral",
        "label": "Neutral",
        "attribute": "neutral",
        "color": ASTROMOL_BLUE,
    },
    {
        "key": "radical",
        "label": "Radical",
        "attribute": "radical",
        "color": "red",
    },
)
MOLECULE_TYPE_BY_SOURCE_MATRIX_TYPE_KEYS = (
    "neutral",
    "radical",
    "cation",
    "anion",
    "cyclic",
)
MOLECULE_TYPE_BY_SOURCE_ENRICHMENT_COLORS = (
    "#6F7377",
    "#F8F8F8",
    ASTROMOL_BLUE,
)
DU_BY_SOURCE_TYPE_SPECS = (
    {
        "key": "carbon_star",
        "label": "Carbon Star",
        "color": "darkorange",
    },
    {
        "key": "dark_cloud",
        "label": "Dark Cloud",
        "color": "forestgreen",
    },
    {
        "key": "diffuse_cloud",
        "label": "Diffuse Cloud",
        "color": "red",
    },
    {
        "key": "sfr",
        "label": "SFR",
        "color": ASTROMOL_BLUE,
    },
)
DU_BY_SOURCE_TYPE_LEGEND_ORDER = (
    "sfr",
    "carbon_star",
    "dark_cloud",
    "diffuse_cloud",
)

FACILITY_SHARE_BACKGROUND = "#F5F6FF"
FACILITY_SHARE_INACTIVE = "#F87070"
FACILITY_SHARE_BAR_SIZE = (7.2, 4.8)
FACILITY_SHARE_BAR_AXES_BOUNDS = (0.23, 0.15, 0.735, 0.82)
LEGACY_SCOPE_COLORS_BY_SHORTNAME = {
    "GBT 100-m": "#377eb8",
    "IRAM 30-m": "#ff7f00",
    "NRAO 140-ft": "#4daf4a",
    "NRAO/ARO 12-m": "#f781bf",
    "NRAO 36-ft": "#a65628",
    "Nobeyama 45-m": "#984ea3",
    "Yebes 40-m": "#999999",
    "ALMA": "#e41a1c",
    "SMT": "#dede00",
}
MODERN_SCOPE_COLORS_BY_LABEL = {
    "NRAO 36-ft": "#B68A72",
    "IRAM 30-m": "#E9A73A",
    "GBT 100-m": "#4E9BCB",
    "ALMA": "#C98CB8",
    "Nobeyama 45-m": "#8E78B7",
    "NRAO/ARO 12-m": "#E78FBC",
    "NRAO 140-ft": "#7FBE7F",
    "SMT": "#BDB765",
    "Yebes 40-m": MIT_RED,
}
MODERN_SCOPE_HIGHLIGHT_LABEL = "Yebes 40-m"
MODERN_SCOPE_DORMANT_LABELS = {
    "NRAO 36-ft",
    "Nobeyama 45-m",
    "NRAO 140-ft",
}
SCOPES_BY_YEAR_CUTOFFS = {
    "NRAO 36-ft": 1985,
    "NRAO 140-ft": 1993,
    "Nobeyama 45-m": 1997,
}
LEGACY_2021_FACILITY_SHARE_NICKS = (
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

ROLLING_RATE_HEATMAP_SIZE = (3.5, 2.8)
ROLLING_RATE_HEATMAP_AXES_BOUNDS = (0.14, 0.145, 0.835, 0.59)
ROLLING_RATE_HEATMAP_COLORBAR_BOUNDS = (0.14, 0.805, 0.835, 0.055)
ROLLING_RATE_HEATMAP_XTICKS = (1940, 1960, 1980, 2000, 2020)
ROLLING_RATE_HEATMAP_LABELS = {
    "2 atoms": "2",
    "3 atoms": "3",
    "4 atoms": "4",
    "5 atoms": "5",
    "6 atoms": "6",
    "7 atoms": "7",
    "8 atoms": "8",
    "9 atoms": "9",
    "10 atoms": "10",
    "11 atoms": "11",
    "12 atoms": "12",
    "13+ atoms": "13+",
    "Fullerenes": "Fuller",
    "PAHs": "PAH",
}

PERIODIC_HEATMAP_ELEMENTS = (
    ("H", 1, 6.9),
    ("He", 18, 6.9),
    ("Li", 1, 5.75),
    ("Be", 2, 5.75),
    ("B", 13, 5.75),
    ("C", 14, 5.75),
    ("N", 15, 5.75),
    ("O", 16, 5.75),
    ("F", 17, 5.75),
    ("Ne", 18, 5.75),
    ("Na", 1, 4.6),
    ("Mg", 2, 4.6),
    ("Al", 13, 4.6),
    ("Si", 14, 4.6),
    ("P", 15, 4.6),
    ("S", 16, 4.6),
    ("Cl", 17, 4.6),
    ("Ar", 18, 4.6),
    ("K", 1, 3.45),
    ("Ca", 2, 3.45),
    ("Sc", 3, 3.45),
    ("Ti", 4, 3.45),
    ("V", 5, 3.45),
    ("Cr", 6, 3.45),
    ("Mn", 7, 3.45),
    ("Fe", 8, 3.45),
    ("Co", 9, 3.45),
    ("Ni", 10, 3.45),
    ("Cu", 11, 3.45),
    ("Zn", 12, 3.45),
    ("Ga", 13, 3.45),
    ("Ge", 14, 3.45),
    ("As", 15, 3.45),
    ("Se", 16, 3.45),
    ("Br", 17, 3.45),
    ("Kr", 18, 3.45),
    ("Rb", 1, 2.3),
    ("Sr", 2, 2.3),
    ("Y", 3, 2.3),
    ("Zr", 4, 2.3),
    ("Nb", 5, 2.3),
    ("Mo", 6, 2.3),
    ("Tc", 7, 2.3),
    ("Ru", 8, 2.3),
    ("Rh", 9, 2.3),
    ("Pd", 10, 2.3),
    ("Ag", 11, 2.3),
    ("Cd", 12, 2.3),
    ("In", 13, 2.3),
    ("Sn", 14, 2.3),
    ("Sb", 15, 2.3),
    ("Te", 16, 2.3),
    ("I", 17, 2.3),
    ("Xe", 18, 2.3),
    ("Cs", 1, 1.15),
    ("Ba", 2, 1.15),
    ("Hf", 4, 1.15),
    ("Ta", 5, 1.15),
    ("W", 6, 1.15),
    ("Re", 7, 1.15),
    ("Os", 8, 1.15),
    ("Ir", 9, 1.15),
    ("Pt", 10, 1.15),
    ("Au", 11, 1.15),
    ("Hg", 12, 1.15),
    ("Tl", 13, 1.15),
    ("Pb", 14, 1.15),
    ("Bi", 15, 1.15),
    ("Po", 16, 1.15),
    ("At", 17, 1.15),
    ("Rn", 18, 1.15),
    ("Fr", 1, 0.0),
    ("Ra", 2, 0.0),
    ("Rf", 4, 0.0),
    ("Db", 5, 0.0),
    ("Sg", 6, 0.0),
    ("Bh", 7, 0.0),
    ("Hs", 8, 0.0),
    ("Mt", 9, 0.0),
    ("Ds", 10, 0.0),
    ("Rg", 11, 0.0),
    ("Cn", 12, 0.0),
    ("Nh", 13, 0.0),
    ("Fl", 14, 0.0),
    ("Mc", 15, 0.0),
    ("Lv", 16, 0.0),
    ("Ts", 17, 0.0),
    ("Og", 18, 0.0),
)
PERIODIC_HEATMAP_START_COLOR = "#f1fb53"
PERIODIC_HEATMAP_STOP_COLOR = "#f00707"
PERIODIC_ELEMENT_FALLBACKS = {
    "Ds": (110, 281.0, "Darmstadtium"),
    "Rg": (111, 280.0, "Roentgenium"),
    "Cn": (112, 285.0, "Copernicium"),
    "Nh": (113, 286.0, "Nihonium"),
    "Fl": (114, 289.0, "Flerovium"),
    "Mc": (115, 290.0, "Moscovium"),
    "Lv": (116, 293.0, "Livermorium"),
    "Ts": (117, 294.0, "Tennessine"),
    "Og": (118, 294.0, "Oganesson"),
}

MASS_BY_WAVELENGTH_SPECS = (
    {
        "key": "cm",
        "label": "centimeter",
        "color": ASTROMOL_BLUE,
        "annotation_xy": (90, 0.01),
        "zorder": 4,
    },
    {
        "key": "mm",
        "label": "millimeter",
        "color": "darkorange",
        "annotation_xy": (60, 0.022),
        "truncate_at_data_max": True,
    },
    {
        "key": "sub-mm",
        "label": "sub-millimeter",
        "color": "forestgreen",
        "annotation_xy": (40, 0.03),
    },
    {
        "key": "IR",
        "label": "infrared",
        "color": "black",
        "annotation_xy": (68, 0.0025),
        "truncate_at_data_max": True,
        "zorder": 5,
    },
    {
        "key": "UV-Vis",
        "label": "visible/ultraviolet",
        "color": "violet",
        "annotation_xy": (14, 0.033),
        "truncate_at_data_max": True,
    },
)
MASS_BY_WAVELENGTH_BOX_SIZE = (7.2, 4.8)
MASS_BY_WAVELENGTH_BOX_COLORS = {
    "cm": ASTROMOL_BLUE,
    "mm": "#E69F00",
    "sub-mm": "#009E73",
    "IR": "black",
    "UV-Vis": "#CC79A7",
}
WAVES_BY_SOURCE_TYPE_WAVELENGTHS = (
    "cm",
    "mm",
    "sub-mm",
    "IR",
    "UV",
    "Vis",
)
WAVES_BY_SOURCE_TYPE_DISPLAY_WAVELENGTHS = {
    "carbon_star": ("cm", "mm", "sub-mm", "IR"),
    "dark_cloud": ("cm", "mm"),
    "diffuse_cloud": ("cm", "mm", "sub-mm", "IR", "UV-Vis"),
    "sfr": ("cm", "mm", "sub-mm"),
}
WAVES_BY_SOURCE_TYPE_PANEL_ORDER = (
    "carbon_star",
    "dark_cloud",
    "diffuse_cloud",
    "sfr",
)
WAVES_BY_SOURCE_TYPE_PANEL_LABELS = {
    "carbon_star": "Carbon Stars",
    "dark_cloud": "Dark Clouds",
    "diffuse_cloud": "Diffuse Clouds",
    "sfr": "SFRs",
}
WAVES_BY_SOURCE_TYPE_COLORS = {
    "cm": ASTROMOL_BLUE,
    "mm": "darkorange",
    "sub-mm": "forestgreen",
    "IR": "black",
    "UV": "violet",
    "Vis": "violet",
    "UV-Vis": "violet",
}
DU_BY_SOURCE_TYPE_BOXPLOT_ORDER = (
    "dark_cloud",
    "carbon_star",
    "sfr",
    "diffuse_cloud",
)
DU_BY_SOURCE_TYPE_BOX_COLORS = {
    "carbon_star": "#E69F00",
    "dark_cloud": "#009E73",
    "sfr": ASTROMOL_BLUE,
    "diffuse_cloud": "#D55E00",
}
RELATIVE_DU_BY_SOURCE_TYPE_PANEL_ORDER = (
    "carbon_star",
    "diffuse_cloud",
    "dark_cloud",
    "sfr",
)
MOLECULES_BY_WAVELENGTH_ATOMS_SPECS = (
    {"key": "cm", "label": "cm", "plot": "kde"},
    {"key": "mm", "label": "mm", "plot": "kde"},
    {"key": "sub-mm", "label": "sub-mm", "plot": "kde"},
    {"key": "IR", "label": "IR", "plot": "kde"},
    {"key": "Vis", "label": "Vis", "plot": "hist"},
    {"key": "UV", "label": "UV", "plot": "hist"},
)
MOLECULES_BY_WAVELENGTH_ATOM_BINS = (
    2,
    3,
    4,
    5,
    6,
    7,
    8,
    9,
    10,
    11,
    12,
    "13+",
)
MOLECULES_BY_WAVELENGTH_BUBBLE_SIZE = (7.2, 4.8)

DETECTION_RATE_BY_ATOMS_CATEGORIES = (
    {"category": 2, "label": "2", "x": 2},
    {"category": 3, "label": "3", "x": 3},
    {"category": 4, "label": "4", "x": 4},
    {"category": 5, "label": "5", "x": 5},
    {"category": 6, "label": "6", "x": 6},
    {"category": 7, "label": "7", "x": 7},
    {"category": 8, "label": "8", "x": 8},
    {"category": 9, "label": "9", "x": 9},
    {"category": 10, "label": "10", "x": 10},
    {"category": 11, "label": "11", "x": 11},
    {"category": 12, "label": "12", "x": 12},
    {"category": "13+", "label": "13+", "x": 13},
    {"category": "PAHs", "label": "PAHs", "x": 15},
    {"category": "Fullerenes", "label": "Fullerenes", "x": 17},
)
DETECTION_RATE_BY_ATOMS_XTICKS = (2, 4, 6, 8, 10, 12, 15, 17)
DETECTION_RATE_BY_ATOMS_XTICK_LABELS = (
    "2",
    "4",
    "6",
    "8",
    "10",
    "12",
    "PAHs",
    "Fullerenes",
)

CUMULATIVE_DETECTION_FACILITY_MARKERS = (
    {
        "label": "NRAO 36-foot (1968)",
        "year": 1968,
        "count_offset": 4,
        "text_offset": (0, 35),
        "rotation": 90,
        "vertical_alignment": "bottom",
        "horizontal_alignment": "center",
    },
    {
        "label": "Nobeyama (1982)",
        "year": 1982,
        "count_offset": -4,
        "text_offset": (0, -35),
        "rotation": 0,
        "vertical_alignment": "top",
        "horizontal_alignment": "left",
    },
    {
        "label": "IRAM (1984)",
        "year": 1984,
        "count_offset": 4,
        "text_offset": (0, 35),
        "rotation": 90,
        "vertical_alignment": "bottom",
        "horizontal_alignment": "center",
    },
    {
        "label": "GBT (2001)",
        "year": 2001,
        "count_offset": 4,
        "text_offset": (0, 35),
        "rotation": 90,
        "vertical_alignment": "bottom",
        "horizontal_alignment": "center",
    },
    {
        "label": "Yebes (2007)",
        "year": 2007,
        "count_offset": -4,
        "text_offset": (0, -35),
        "rotation": 90,
        "vertical_alignment": "top",
        "horizontal_alignment": "center",
    },
    {
        "label": "ALMA (2011)",
        "year": 2011,
        "count_offset": -4,
        "text_offset": (0, -35),
        "rotation": 90,
        "vertical_alignment": "top",
        "horizontal_alignment": "center",
    },
)

LEGACY_CUMULATIVE_BY_ATOMS_SERIES = (
    {"category": 2, "label": "2 atoms", "color": "#000000"},
    {"category": 3, "label": "3 atoms", "color": "#800000"},
    {"category": 4, "label": "4 atoms", "color": "#f032e6"},
    {"category": 5, "label": "5 atoms", "color": "#9A6324"},
    {"category": 6, "label": "6 atoms", "color": "dodgerblue"},
    {"category": 7, "label": "7 atoms", "color": "#e6194B"},
    {"category": 8, "label": "8 atoms", "color": "#469990"},
    {"category": 9, "label": "9 atoms", "color": "#f58231"},
    {"category": 10, "label": "10 atoms", "color": "#42d4f4"},
    {"category": 11, "label": "11 atoms", "color": "#ffe119"},
    {"category": 12, "label": "12 atoms", "color": "#3cb44b"},
    {"category": "13+", "label": "13+ atoms", "color": "#e6beff"},
    {"category": "Fullerenes", "label": "Fullerenes", "color": "#000075"},
    {"category": "PAHs", "label": "PAHs", "color": "#aaffc3"},
)

COLORBLIND_CUMULATIVE_BY_ATOMS_SERIES = (
    {"category": 2, "label": "2 atoms", "color": "#000000"},
    {"category": 3, "label": "3 atoms", "color": "#D55E00"},
    {"category": 4, "label": "4 atoms", "color": "#0072B2"},
    {"category": 5, "label": "5 atoms", "color": "#009E73"},
    {"category": 6, "label": "6 atoms", "color": "#CC79A7"},
    {"category": 7, "label": "7 atoms", "color": "#E69F00"},
    {"category": 8, "label": "8 atoms", "color": "#56B4E9"},
    {"category": 9, "label": "9 atoms", "color": "#7F7F7F"},
    {"category": 10, "label": "10 atoms", "color": "#882255"},
    {"category": 11, "label": "11 atoms", "color": "#44AA99"},
    {"category": 12, "label": "12 atoms", "color": "#999933"},
    {"category": "13+", "label": "13+ atoms", "color": "#AA4499"},
    {"category": "Fullerenes", "label": "Fullerenes", "color": "#332288"},
    {"category": "PAHs", "label": "PAHs", "color": "#117733"},
)

CUMULATIVE_BY_ATOMS_SERIES = COLORBLIND_CUMULATIVE_BY_ATOMS_SERIES


@dataclass(frozen=True)
class DetectionTrend:
    """Linear trend fitted to a section of a cumulative detection curve."""

    label: str
    start_year: int
    stop_year: int | None
    slope: float


@dataclass(frozen=True)
class CumulativeDetectionData:
    """Data needed to plot cumulative detections over time."""

    years: np.ndarray
    counts: np.ndarray
    first_detection_years: dict[str, int]
    trends: tuple[DetectionTrend, ...]
    show_facility_markers: bool = False

    @property
    def total(self) -> int:
        """Final cumulative detection count."""
        return int(self.counts[-1])

    @property
    def start_year(self) -> int:
        """First year in the cumulative series."""
        return int(self.years[0])

    @property
    def end_year(self) -> int:
        """Last year in the cumulative series."""
        return int(self.years[-1])


@dataclass(frozen=True)
class CumulativeByAtomsSeries:
    """One cumulative-detection series for an atom-count category."""

    category: int | str
    label: str
    color: str
    years: np.ndarray
    counts: np.ndarray
    first_detection_years: dict[str, int]

    @property
    def final_count(self) -> int:
        """Final count for this series."""
        return int(self.counts[-1])


@dataclass(frozen=True)
class CumulativeByAtomsData:
    """Data needed to plot cumulative detections grouped by atom count."""

    years: np.ndarray
    series: tuple[CumulativeByAtomsSeries, ...]

    @property
    def start_year(self) -> int:
        """First year in the cumulative series."""
        return int(self.years[0])

    @property
    def end_year(self) -> int:
        """Last year in the cumulative series."""
        return int(self.years[-1])

    @property
    def max_count(self) -> int:
        """Largest cumulative count across all plotted series."""
        return max(series.final_count for series in self.series)

    @property
    def total(self) -> int:
        """Total number of molecules counted across all series."""
        return sum(series.final_count for series in self.series)


@dataclass(frozen=True)
class RollingRateHeatmapData:
    """Data needed to plot rolling detection rates by atom-count category."""

    years: np.ndarray
    labels: tuple[str, ...]
    matrix: np.ndarray
    window: int
    vmax: float

    @property
    def start_year(self) -> int:
        """First year in the rolling-rate series."""
        return int(self.years[0])

    @property
    def end_year(self) -> int:
        """Last year in the rolling-rate series."""
        return int(self.years[-1])


@dataclass(frozen=True)
class PeriodicHeatmapCell:
    """One displayed cell in the periodic-table heatmap."""

    symbol: str
    atomic_number: int
    atomic_mass: float
    name: str
    group: int
    y_position: float
    count: int

    @property
    def x_position(self) -> float:
        """Left x-coordinate for this periodic-table cell."""
        return float(self.group - 1)


@dataclass(frozen=True)
class PeriodicHeatmapData:
    """Data needed to plot the elemental-composition periodic heatmap."""

    cells: tuple[PeriodicHeatmapCell, ...]
    element_counts: dict[str, int]
    molecule_count: int

    @property
    def detected_element_count(self) -> int:
        """Number of elements present in one or more molecules."""
        return len(self.element_counts)

    @property
    def max_count(self) -> int:
        """Largest molecule count for any element."""
        return max(self.element_counts.values())


@dataclass(frozen=True)
class DUHistogramData:
    """Data needed to plot degree-of-unsaturation histograms."""

    values: tuple[float, ...]
    molecule_labels: tuple[str, ...]
    formula_labels: tuple[str, ...]
    molecule_count: int

    @property
    def saturated_count(self) -> int:
        """Number of included molecules with DU = 0."""
        return sum(value == 0 for value in self.values)

    @property
    def unsaturated_count(self) -> int:
        """Number of included molecules with DU > 0."""
        return sum(value > 0 for value in self.values)

    @property
    def unsaturated_percent(self) -> int:
        """Integer unsaturated percentage, matching legacy truncation."""
        if not self.values:
            return 0
        return int(100 * self.unsaturated_count / len(self.values))

    @property
    def max_du(self) -> float:
        """Largest degree of unsaturation in the sample."""
        return max(self.values)

    def histogram_counts(
        self,
        bins: tuple[float, ...] = DU_HISTOGRAM_BINS,
    ) -> np.ndarray:
        """Return counts in the legacy half-DU-width bins."""
        counts, _ = np.histogram(self.values, bins=np.array(bins))
        return counts

    @property
    def counts_by_value(self) -> dict[float, int]:
        """Return exact DU-value counts."""
        return self.value_counts()

    def value_counts(self, *, include_negative_du: bool = True) -> dict[float, int]:
        """Return exact DU-value counts, optionally excluding domain artifacts."""
        counts = defaultdict(int)
        for value in self.values:
            if value < 0 and not include_negative_du:
                continue
            counts[value] += 1
        return dict(counts)


@dataclass(frozen=True)
class KappaHistogramData:
    """Data needed to plot Ray asymmetry-parameter histograms."""

    values: tuple[float, ...]
    molecule_labels: tuple[str, ...]
    molecule_count: int

    @property
    def min_kappa(self) -> float:
        """Smallest Ray asymmetry parameter in the sample."""
        return min(self.values)

    @property
    def max_kappa(self) -> float:
        """Largest Ray asymmetry parameter in the sample."""
        return max(self.values)

    def histogram_counts(self, bins: int = 100) -> np.ndarray:
        """Return histogram counts using the legacy integer-bin convention."""
        counts, _ = np.histogram(self.values, bins=bins)
        return counts


@dataclass(frozen=True)
class MoleculeTypeCategory:
    """One molecule-type category used by the type ring chart."""

    key: str
    label: str
    plural_label: str
    color: str
    percent_color: str
    count: int
    fraction: float

    @property
    def percent(self) -> float:
        """Category percentage of the molecule sample."""
        return 100.0 * self.fraction


@dataclass(frozen=True)
class MoleculeTypeData:
    """Data needed to plot the molecule-type ring chart."""

    categories: tuple[MoleculeTypeCategory, ...]
    molecule_count: int

    @property
    def counts(self) -> dict[str, int]:
        """Category counts keyed by short category name."""
        return {category.key: category.count for category in self.categories}

    @property
    def fractions(self) -> dict[str, float]:
        """Category fractions keyed by short category name."""
        return {category.key: category.fraction for category in self.categories}


@dataclass(frozen=True)
class SourceTypeCategory:
    """One generalized first-detection source category."""

    key: str
    label: str
    color: str
    count: int
    fraction: float
    label_y: float
    percent_y: float

    @property
    def percent(self) -> float:
        """Category percentage of the molecule sample."""
        return 100.0 * self.fraction


@dataclass(frozen=True)
class SourceTypeData:
    """Data needed to plot the generalized source-type ring chart."""

    categories: tuple[SourceTypeCategory, ...]
    molecule_count: int

    @property
    def counts(self) -> dict[str, int]:
        """Generalized source-type counts keyed by category name."""
        return {category.key: category.count for category in self.categories}

    @property
    def fractions(self) -> dict[str, float]:
        """Generalized source-type fractions keyed by category name."""
        return {category.key: category.fraction for category in self.categories}


@dataclass(frozen=True)
class IndividualSourceData:
    """Data needed to plot the individual first-detection source ring chart."""

    categories: tuple[SourceTypeCategory, ...]
    molecule_count: int

    @property
    def counts(self) -> dict[str, int]:
        """Individual-source contribution counts keyed by category name."""
        return {category.key: category.count for category in self.categories}

    @property
    def fractions(self) -> dict[str, float]:
        """Individual-source contribution fractions keyed by category name."""
        return {category.key: category.fraction for category in self.categories}


@dataclass(frozen=True)
class MoleculeTypesForSource:
    """Molecule-type counts credited to one generalized first-detection source."""

    key: str
    label: str
    source_count: int
    counts_by_type: dict[str, int]

    @property
    def total_type_count(self) -> int:
        """Total molecule-type credits in this source category."""
        return sum(self.counts_by_type.values())


@dataclass(frozen=True)
class MoleculeTypeBySourceData:
    """Data needed to plot molecule-type mixes by source category."""

    categories: tuple[MoleculeTypesForSource, ...]
    molecule_count: int
    overall_type_counts: dict[str, int]

    @property
    def counts(self) -> dict[str, dict[str, int]]:
        """Nested counts keyed first by source category and then molecule type."""
        return {
            category.key: dict(category.counts_by_type)
            for category in self.categories
        }

    @property
    def source_counts(self) -> dict[str, int]:
        """Source-category molecule counts keyed by source category."""
        return {
            category.key: category.source_count
            for category in self.categories
        }


@dataclass(frozen=True)
class DUBySourceTypeCategory:
    """DU values credited to one generalized first-detection source category."""

    key: str
    label: str
    color: str
    values: tuple[float, ...]

    @property
    def count(self) -> int:
        """Number of DU values in this source category."""
        return len(self.values)


@dataclass(frozen=True)
class DUBySourceTypeData:
    """Data needed to plot DU distributions by generalized source category."""

    categories: tuple[DUBySourceTypeCategory, ...]
    molecule_count: int

    @property
    def counts(self) -> dict[str, int]:
        """DU-value counts keyed by generalized source category."""
        return {
            category.key: category.count
            for category in self.categories
        }

    @property
    def values_by_source(self) -> dict[str, tuple[float, ...]]:
        """DU values keyed by generalized source category."""
        return {
            category.key: category.values
            for category in self.categories
        }


@dataclass(frozen=True)
class _DUBySourceTypeCurve:
    """Rendered DU/source KDE curve data used for count-label placement."""

    category: DUBySourceTypeCategory
    x_values: np.ndarray
    density: np.ndarray

    @property
    def peak_x(self) -> float:
        """x-coordinate of the highest KDE value."""
        return float(self.x_values[int(np.argmax(self.density))])

    @property
    def peak_y(self) -> float:
        """Highest KDE value."""
        return float(np.max(self.density))


@dataclass(frozen=True)
class RelativeDUBySourceTypeCategory:
    """Relative-DU values credited to one generalized source category."""

    key: str
    label: str
    color: str
    values: tuple[float, ...]

    @property
    def count(self) -> int:
        """Number of relative-DU values in this source category."""
        return len(self.values)


@dataclass(frozen=True)
class RelativeDUBySourceTypeData:
    """Data needed to plot relative-DU distributions by source category."""

    categories: tuple[RelativeDUBySourceTypeCategory, ...]
    molecule_count: int

    @property
    def counts(self) -> dict[str, int]:
        """Relative-DU value counts keyed by source category."""
        return {
            category.key: category.count
            for category in self.categories
        }

    @property
    def values_by_source(self) -> dict[str, tuple[float, ...]]:
        """Relative-DU values keyed by source category."""
        return {
            category.key: category.values
            for category in self.categories
        }


@dataclass(frozen=True)
class MassBySourceTypeCategory:
    """Molecular masses credited to one generalized source category."""

    key: str
    label: str
    color: str
    masses: tuple[float, ...]

    @property
    def count(self) -> int:
        """Number of molecular masses in this source category."""
        return len(self.masses)


@dataclass(frozen=True)
class MassBySourceTypeData:
    """Data needed to plot molecular-mass distributions by source category."""

    categories: tuple[MassBySourceTypeCategory, ...]
    molecule_count: int
    mass_range: tuple[float, float]

    @property
    def counts(self) -> dict[str, int]:
        """Molecular-mass counts keyed by generalized source category."""
        return {
            category.key: category.count
            for category in self.categories
        }

    @property
    def masses_by_source(self) -> dict[str, tuple[float, ...]]:
        """Molecular masses keyed by generalized source category."""
        return {
            category.key: category.masses
            for category in self.categories
        }


@dataclass(frozen=True)
class WavelengthBySourceTypeCategory:
    """First-detection wavelength counts for one source category."""

    key: str
    label: str
    counts_by_wavelength: dict[str, int]

    @property
    def count(self) -> int:
        """Total wavelength credits in this source category."""
        return sum(self.counts_by_wavelength.values())

    def count_for_display_wavelength(self, wavelength: str) -> int:
        """Return counts for a displayed wavelength, combining UV/Vis if needed."""
        if wavelength == "UV-Vis":
            return (
                self.counts_by_wavelength.get("UV", 0)
                + self.counts_by_wavelength.get("Vis", 0)
            )
        return self.counts_by_wavelength.get(wavelength, 0)


@dataclass(frozen=True)
class WavelengthBySourceTypeData:
    """Data needed to plot wavelength mixes by source category."""

    categories: tuple[WavelengthBySourceTypeCategory, ...]
    molecule_count: int

    @property
    def counts(self) -> dict[str, dict[str, int]]:
        """Nested counts keyed first by source category, then wavelength."""
        return {
            category.key: dict(category.counts_by_wavelength)
            for category in self.categories
        }


def _du_histogram_bins(data: DUHistogramData) -> tuple[float, ...]:
    """Return half-DU-width bins covering the legacy range or the data range."""
    max_du = max(12.0, np.ceil(data.max_du * 2) / 2)
    return tuple(np.arange(-0.25, max_du + 0.5, 0.5))


def _matplotlib_formula_text(formula: str) -> str:
    """Return a compact formula label with simple mathtext subscripts."""
    import re

    return re.sub(r"(?<=[A-Za-z\]])(\d+)", r"$_{\1}$", formula)


def _matplotlib_formula_list(formulas: list[str]) -> str:
    """Return compact display text for one or more formula labels."""
    unique_formulas = list(dict.fromkeys(formulas))
    split_formulas = [
        formula.split("-", 1)
        for formula in unique_formulas
        if "-" in formula
    ]
    if len(split_formulas) == len(unique_formulas):
        suffixes = {suffix for _, suffix in split_formulas}
        if len(suffixes) == 1:
            suffix = split_formulas[0][1]
            prefixes = "/".join(f"{prefix}-" for prefix, _ in split_formulas)
            return f"{prefixes}{_matplotlib_formula_text(suffix)}"
    return ", ".join(
        _matplotlib_formula_text(formula)
        for formula in unique_formulas
    )


@dataclass(frozen=True)
class MassByWavelengthSeries:
    """Molecular masses contributing to one wavelength KDE trace."""

    key: str
    label: str
    color: str
    masses: tuple[float, ...]
    annotation_xy: tuple[float, float]
    truncate_at_data_max: bool = False
    zorder: int | None = None

    @property
    def count(self) -> int:
        """Number of molecule masses in the wavelength sample."""
        return len(self.masses)


@dataclass(frozen=True)
class MassByWavelengthData:
    """Data needed to plot molecular mass distributions by wavelength."""

    series: tuple[MassByWavelengthSeries, ...]
    molecule_count: int

    @property
    def counts(self) -> dict[str, int]:
        """Sample counts keyed by wavelength label."""
        return {series.key: series.count for series in self.series}


@dataclass(frozen=True)
class _MassByWavelengthCurve:
    """Rendered KDE curve data used for automatic count-label placement."""

    series: MassByWavelengthSeries
    x_values: np.ndarray
    density: np.ndarray

    @property
    def peak_x(self) -> float:
        """x-coordinate of the highest KDE value."""
        return float(self.x_values[int(np.argmax(self.density))])

    @property
    def peak_y(self) -> float:
        """Highest KDE value."""
        return float(np.max(self.density))


@dataclass(frozen=True)
class MoleculesByWavelengthAtomsSeries:
    """Atom counts contributing to one wavelength-panel distribution."""

    key: str
    label: str
    plot: str
    atom_counts: tuple[int, ...]

    @property
    def count(self) -> int:
        """Number of molecules in this wavelength sample."""
        return len(self.atom_counts)


@dataclass(frozen=True)
class MoleculesByWavelengthAtomsData:
    """Data needed to plot atom-count distributions by detection wavelength."""

    series: tuple[MoleculesByWavelengthAtomsSeries, ...]
    molecule_count: int

    @property
    def counts(self) -> dict[str, int]:
        """Sample counts keyed by wavelength label."""
        return {series.key: series.count for series in self.series}

    @property
    def max_atoms(self) -> int:
        """Largest atom count in any wavelength sample."""
        return max(
            atom_count
            for series in self.series
            for atom_count in series.atom_counts
        )

    def matrix(
        self,
        bins: tuple[int | str, ...] = MOLECULES_BY_WAVELENGTH_ATOM_BINS,
    ) -> np.ndarray:
        """Return wavelength-by-atom-count matrix for discrete heatmaps."""
        values = np.zeros((len(self.series), len(bins)), dtype=int)
        for row, series in enumerate(self.series):
            for atom_count in series.atom_counts:
                if isinstance(bins[-1], str) and atom_count >= int(bins[-2]) + 1:
                    values[row, len(bins) - 1] += 1
                elif atom_count in bins:
                    values[row, bins.index(atom_count)] += 1
        return values


@dataclass(frozen=True)
class DetectionRateByAtomsPoint:
    """One average detection-rate point for an atom-count category."""

    category: int | str
    label: str
    x_position: float
    count: int
    first_year: int | None
    rate: float


@dataclass(frozen=True)
class DetectionRateByAtomsData:
    """Data needed to plot average detection rate by atom-count category."""

    end_year: int
    points: tuple[DetectionRateByAtomsPoint, ...]

    @property
    def visible_points(self) -> tuple[DetectionRateByAtomsPoint, ...]:
        """Return points visible in the legacy 2-12 atom x-axis range."""
        return tuple(point for point in self.points if 1 <= point.x_position <= 12.8)


@dataclass(frozen=True)
class FacilityShare:
    """One facility's share of first detections during its operating window."""

    nick: str
    label: str
    start_year: int
    end_year: int
    detection_count: int
    total_window_detections: int
    active_at_view: bool

    @property
    def fraction(self) -> float:
        """Facility contribution fraction during its operating window."""
        if self.total_window_detections == 0:
            return 0.0
        return self.detection_count / self.total_window_detections

    @property
    def percent(self) -> int:
        """Integer percentage label, matching legacy truncation."""
        return int(self.fraction * 100)


@dataclass(frozen=True)
class FacilityShareData:
    """Data needed to plot facility shares."""

    end_year: int
    facilities: tuple[FacilityShare, ...]


@dataclass(frozen=True)
class ScopeDetectionSeries:
    """Cumulative first-detection contribution series for one facility."""

    nick: str
    label: str
    start_year: int
    fit_stop_year: int
    color: str
    years: np.ndarray
    counts: np.ndarray
    detection_count: int
    rate: float

    @property
    def final_count(self) -> int:
        """Final cumulative contribution count."""
        return int(self.counts[-1])


@dataclass(frozen=True)
class ScopesByYearData:
    """Data needed to plot cumulative facility contributions over time."""

    years: np.ndarray
    series: tuple[ScopeDetectionSeries, ...]
    end_year: int

    @property
    def max_count(self) -> int:
        """Largest final cumulative contribution count."""
        return max(series.final_count for series in self.series)


def first_detections_by_molecule(
    view: CensusView,
    detection_type: str = "ISM/CSM",
    *,
    include_tentative: bool = False,
    include_disputed: bool = False,
    include_isotopologues: bool = False,
) -> dict[str, Detection]:
    """Return the first detection record for each molecule in a view."""
    detections_by_molecule = defaultdict(list)
    for detection in view.context_detections(
        detection_type,
        include_tentative=include_tentative,
        include_disputed=include_disputed,
        include_isotopologues=include_isotopologues,
    ):
        detections_by_molecule[detection.molecule.label].append(detection)

    return {
        label: min(
            detections,
            key=lambda detection: (
                detection.sortdate,
                detection.id,
            ),
        )
        for label, detections in detections_by_molecule.items()
    }


def first_detection_years_by_molecule(
    view: CensusView,
    detection_type: str = "ISM/CSM",
    *,
    include_tentative: bool = False,
    include_disputed: bool = False,
    include_isotopologues: bool = False,
) -> dict[str, int]:
    """Return the first detection year for each molecule in a view."""
    return {
        label: detection.year
        for label, detection in first_detections_by_molecule(
            view,
            detection_type,
            include_tentative=include_tentative,
            include_disputed=include_disputed,
            include_isotopologues=include_isotopologues,
        ).items()
    }


def _detection_matches_atom_category(
    detection: Detection,
    category: int | str,
) -> bool:
    """Return whether a detection belongs to an atom-count figure category."""
    molecule = detection.molecule
    if isinstance(category, int):
        return molecule.natoms == category
    if category == "13+":
        return (
            molecule.natoms >= 13
            and not molecule.pah
            and not molecule.fullerene
        )
    if category == "PAHs":
        return molecule.pah
    if category == "Fullerenes":
        return molecule.fullerene
    raise ValueError(f"Unknown atom-count category: {category}")


def _counts_for_detection_years(
    detection_years: list[int],
    years: np.ndarray,
) -> np.ndarray:
    """Return cumulative counts for detection years over fixed plot years."""
    return np.array(
        [
            sum(detection_year <= year for detection_year in detection_years)
            for year in years
        ]
    )


def cumulative_by_atoms_series_specs(
    view: CensusView,
) -> tuple[dict[str, object], ...]:
    """Return the default atom-count palette for a census/current view.

    Legacy 2018/2021 reproduction views keep the original census palette.
    Current and future-facing views use the color-blind-friendly palette.
    """
    if not view.is_current and view.census_year is not None and view.census_year <= 2021:
        return LEGACY_CUMULATIVE_BY_ATOMS_SERIES
    return COLORBLIND_CUMULATIVE_BY_ATOMS_SERIES


def cumulative_by_atoms_data(
    view: CensusView,
    detection_type: str = "ISM/CSM",
    *,
    start_year: int | None = None,
    end_year: int | None = None,
    include_tentative: bool = False,
    include_disputed: bool = False,
    include_isotopologues: bool = False,
    series_specs: tuple[dict[str, object], ...] | None = None,
) -> CumulativeByAtomsData:
    """Return cumulative-detection data grouped by atom-count category.

    The default grouping follows the legacy figure: 2-12 atoms, a true 13+
    category excluding PAHs and fullerenes, and separate PAH/fullerene traces.
    """
    first_detections = first_detections_by_molecule(
        view,
        detection_type,
        include_tentative=include_tentative,
        include_disputed=include_disputed,
        include_isotopologues=include_isotopologues,
    )
    if not first_detections:
        raise ValueError("Cannot build an atom-count plot with no detections.")
    if series_specs is None:
        series_specs = cumulative_by_atoms_series_specs(view)

    if start_year is None:
        start_year = min(detection.year for detection in first_detections.values())
    if end_year is None:
        end_year = view.census_year if not view.is_current else date.today().year
    if end_year < start_year:
        raise ValueError("end_year must be greater than or equal to start_year.")

    years = np.arange(start_year, end_year + 1)
    series = []
    for spec in series_specs:
        category = spec["category"]
        category_detections = {
            label: detection
            for label, detection in first_detections.items()
            if _detection_matches_atom_category(detection, category)
        }
        first_years = {
            label: detection.year
            for label, detection in category_detections.items()
        }
        series.append(
            CumulativeByAtomsSeries(
                category=category,
                label=spec["label"],
                color=spec["color"],
                years=years,
                counts=_counts_for_detection_years(
                    list(first_years.values()),
                    years,
                ),
                first_detection_years=first_years,
            )
        )

    return CumulativeByAtomsData(years=years, series=tuple(series))


def annual_detection_counts(
    first_detection_years: dict[str, int],
    years: np.ndarray,
) -> np.ndarray:
    """Return annual first-detection counts over fixed plot years."""
    counts_by_year = defaultdict(int)
    for detection_year in first_detection_years.values():
        counts_by_year[detection_year] += 1

    return np.array([counts_by_year[int(year)] for year in years])


def trailing_rolling_detection_rate(
    first_detection_years: dict[str, int],
    years: np.ndarray,
    *,
    window: int = 10,
) -> np.ndarray:
    """Return trailing rolling detection rates in detections per year.

    A 10-year value at 2026 counts detections from 2017 through 2026 and
    divides by 10. The fixed denominator keeps rates comparable across years.
    """
    if window < 1:
        raise ValueError("window must be a positive integer.")

    annual_counts = annual_detection_counts(first_detection_years, years)
    rates = []
    for index in range(len(years)):
        window_start = max(0, index - window + 1)
        rates.append(annual_counts[window_start : index + 1].sum() / window)
    return np.array(rates)


def rolling_rate_by_atoms_heatmap_data(
    view: CensusView,
    detection_type: str = "ISM/CSM",
    *,
    start_year: int | None = None,
    end_year: int | None = None,
    include_tentative: bool = False,
    include_disputed: bool = False,
    include_isotopologues: bool = False,
    series_specs: tuple[dict[str, object], ...] | None = None,
    window: int = 10,
    vmax: float = 2.0,
) -> RollingRateHeatmapData:
    """Return rolling detection-rate heatmap data grouped by atom count."""
    cumulative_data = cumulative_by_atoms_data(
        view,
        detection_type,
        start_year=start_year,
        end_year=end_year,
        include_tentative=include_tentative,
        include_disputed=include_disputed,
        include_isotopologues=include_isotopologues,
        series_specs=series_specs,
    )

    labels = tuple(
        ROLLING_RATE_HEATMAP_LABELS[series.label]
        for series in cumulative_data.series
    )
    matrix = np.array(
        [
            trailing_rolling_detection_rate(
                series.first_detection_years,
                cumulative_data.years,
                window=window,
            )
            for series in cumulative_data.series
        ]
    )
    return RollingRateHeatmapData(
        years=cumulative_data.years,
        labels=labels,
        matrix=matrix,
        window=window,
        vmax=vmax,
    )


def periodic_heatmap_data(
    view: CensusView,
    detection_type: str = "ISM/CSM",
    *,
    include_tentative: bool = False,
    include_disputed: bool = False,
    include_isotopologues: bool = False,
) -> PeriodicHeatmapData:
    """Return element-composition counts for the periodic-table heatmap."""
    molecules = view.context_molecules(
        detection_type,
        include_tentative=include_tentative,
        include_disputed=include_disputed,
        include_isotopologues=include_isotopologues,
    )
    if not molecules:
        raise ValueError("Cannot build a periodic heatmap with no molecules.")

    element_counts = defaultdict(int)
    for molecule in molecules:
        for symbol, count in molecule.atom_counts.items():
            if count > 0:
                element_counts[symbol] += 1

    cells = []
    for symbol, group, y_position in PERIODIC_HEATMAP_ELEMENTS:
        try:
            element = ELEMENTS[symbol]
            atomic_number = element.number
            atomic_mass = element.mass
            name = element.name.capitalize()
        except KeyError:
            atomic_number, atomic_mass, name = PERIODIC_ELEMENT_FALLBACKS[symbol]
        cells.append(
            PeriodicHeatmapCell(
                symbol=symbol,
                atomic_number=atomic_number,
                atomic_mass=atomic_mass,
                name=name,
                group=group,
                y_position=y_position,
                count=element_counts.get(symbol, 0),
            )
        )

    return PeriodicHeatmapData(
        cells=tuple(cells),
        element_counts=dict(element_counts),
        molecule_count=len(molecules),
    )


def du_histogram_data(
    view: CensusView,
    detection_type: str = "ISM/CSM",
    *,
    include_tentative: bool = False,
    include_disputed: bool = False,
    include_isotopologues: bool = False,
    include_fullerenes: bool = False,
) -> DUHistogramData:
    """Return degree-of-unsaturation values for the legacy histogram domain."""
    molecules = view.context_molecules(
        detection_type,
        include_tentative=include_tentative,
        include_disputed=include_disputed,
        include_isotopologues=include_isotopologues,
    )
    if not molecules:
        raise ValueError("Cannot build a DU histogram with no molecules.")

    values = []
    molecule_labels = []
    formula_labels = []
    for molecule in molecules:
        if molecule.fullerene and not include_fullerenes:
            continue
        if not set(molecule.atoms).issubset(DU_HISTOGRAM_ALLOWED_ELEMENTS):
            continue
        if molecule.du is None:
            continue
        values.append(float(molecule.du))
        molecule_labels.append(molecule.label)
        formula_labels.append(molecule.table_formula)

    if not values:
        raise ValueError("No molecules in this view have DU-compatible formulas.")

    return DUHistogramData(
        values=tuple(values),
        molecule_labels=tuple(molecule_labels),
        formula_labels=tuple(formula_labels),
        molecule_count=len(values),
    )


def kappa_histogram_data(
    view: CensusView,
    detection_type: str = "ISM/CSM",
    *,
    include_tentative: bool = False,
    include_disputed: bool = False,
    include_isotopologues: bool = False,
) -> KappaHistogramData:
    """Return Ray asymmetry-parameter values for molecules with A/B/C constants."""
    molecules = view.context_molecules(
        detection_type,
        include_tentative=include_tentative,
        include_disputed=include_disputed,
        include_isotopologues=include_isotopologues,
    )
    if not molecules:
        raise ValueError("Cannot build a kappa histogram with no molecules.")

    values = []
    molecule_labels = []
    for molecule in molecules:
        if molecule.kappa is None:
            continue
        values.append(float(molecule.kappa))
        molecule_labels.append(molecule.label)

    if not values:
        raise ValueError("No molecules in this view have Ray kappa values.")

    return KappaHistogramData(
        values=tuple(values),
        molecule_labels=tuple(molecule_labels),
        molecule_count=len(values),
    )


def molecule_type_data(
    view: CensusView,
    detection_type: str = "ISM/CSM",
    *,
    include_tentative: bool = False,
    include_disputed: bool = False,
    include_isotopologues: bool = False,
) -> MoleculeTypeData:
    """Return molecule-type category counts for a census view.

    Categories are counted independently, so one molecule may contribute to
    more than one category. This preserves the legacy manuscript convention
    used for neutral radical species, cyclic ions, PAHs, and fullerenes.
    """
    molecules = view.context_molecules(
        detection_type,
        include_tentative=include_tentative,
        include_disputed=include_disputed,
        include_isotopologues=include_isotopologues,
    )
    if not molecules:
        raise ValueError("Cannot build a molecule-type chart with no molecules.")

    molecule_count = len(molecules)
    categories = []
    for spec in TYPE_PIE_SPECS:
        count = sum(
            getattr(molecule, spec["attribute"]) is True
            for molecule in molecules
        )
        categories.append(
            MoleculeTypeCategory(
                key=spec["key"],
                label=spec["label"],
                plural_label=spec["plural_label"],
                color=spec["color"],
                percent_color=spec["percent_color"],
                count=count,
                fraction=count / molecule_count,
            )
        )

    return MoleculeTypeData(
        categories=tuple(categories),
        molecule_count=molecule_count,
    )


def _source_type_category(source_type: str) -> str:
    """Return the generalized source-pie category for one source type."""
    if source_type in SOURCE_PIE_DIRECT_SOURCE_TYPES:
        for spec in SOURCE_PIE_SPECS:
            if source_type in spec["source_types"]:
                return spec["key"]
    return "other"


def source_type_data(
    view: CensusView,
    detection_type: str = "ISM/CSM",
    *,
    include_tentative: bool = False,
    include_disputed: bool = False,
    include_isotopologues: bool = False,
) -> SourceTypeData:
    """Return generalized first-detection source-type counts for a view.

    Each molecule can receive credit for multiple source categories if its
    first detection lists sources in multiple categories, but it is credited at
    most once per generalized category.
    """
    first_detections = first_detections_by_molecule(
        view,
        detection_type,
        include_tentative=include_tentative,
        include_disputed=include_disputed,
        include_isotopologues=include_isotopologues,
    )
    if not first_detections:
        raise ValueError("Cannot build a source-type chart with no detections.")

    counts = defaultdict(int)
    for detection in first_detections.values():
        credited_categories = {
            _source_type_category(source.type)
            for source in detection.sources
        }
        for category_key in credited_categories:
            counts[category_key] += 1

    molecule_count = len(first_detections)
    categories = []
    for spec in SOURCE_PIE_SPECS:
        count = counts[spec["key"]]
        categories.append(
            SourceTypeCategory(
                key=spec["key"],
                label=spec["label"],
                color=spec["color"],
                count=count,
                fraction=count / molecule_count,
                label_y=spec["label_y"],
                percent_y=spec["percent_y"],
            )
        )

    return SourceTypeData(
        categories=tuple(categories),
        molecule_count=molecule_count,
    )


def individual_source_data(
    view: CensusView,
    detection_type: str = "ISM/CSM",
    *,
    include_tentative: bool = False,
    include_disputed: bool = False,
    include_isotopologues: bool = False,
) -> IndividualSourceData:
    """Return individual first-detection source contribution counts.

    This preserves the legacy convention: each listed source on a first
    detection contributes one count. The four named sources are counted
    individually; every other source contribution is counted in ``Other``.
    The denominator remains the number of first-detected molecules, so
    percentages may sum to more than 100%.
    """
    first_detections = first_detections_by_molecule(
        view,
        detection_type,
        include_tentative=include_tentative,
        include_disputed=include_disputed,
        include_isotopologues=include_isotopologues,
    )
    if not first_detections:
        raise ValueError("Cannot build an individual-source chart with no detections.")

    counts = defaultdict(int)
    for detection in first_detections.values():
        for source in detection.sources:
            category_key = INDIVIDUAL_SOURCE_NICK_TO_KEY.get(source.nick, "other")
            counts[category_key] += 1

    molecule_count = len(first_detections)
    categories = []
    for spec in INDIVIDUAL_SOURCE_PIE_SPECS:
        count = counts[spec["key"]]
        categories.append(
            SourceTypeCategory(
                key=spec["key"],
                label=spec["label"],
                color=spec["color"],
                count=count,
                fraction=count / molecule_count,
                label_y=0.0,
                percent_y=0.0,
            )
        )

    return IndividualSourceData(
        categories=tuple(categories),
        molecule_count=molecule_count,
    )


def molecule_type_by_source_type_data(
    view: CensusView,
    detection_type: str = "ISM/CSM",
    *,
    include_tentative: bool = False,
    include_disputed: bool = False,
    include_isotopologues: bool = False,
) -> MoleculeTypeBySourceData:
    """Return molecule-type counts for generalized first-detection sources.

    Source categories are credited once per molecule. Molecule-type categories
    are counted independently within each source category, preserving the
    legacy convention that a molecule can contribute to multiple type wedges.
    """
    first_detections = first_detections_by_molecule(
        view,
        detection_type,
        include_tentative=include_tentative,
        include_disputed=include_disputed,
        include_isotopologues=include_isotopologues,
    )
    if not first_detections:
        raise ValueError(
            "Cannot build a molecule-type-by-source chart with no detections."
        )

    source_keys = {
        spec["key"]
        for spec in MOLECULE_TYPE_BY_SOURCE_SOURCE_SPECS
    }
    counts_by_source = {
        source_key: defaultdict(int)
        for source_key in source_keys
    }
    source_counts = defaultdict(int)
    overall_type_counts = defaultdict(int)

    for detection in first_detections.values():
        for type_spec in MOLECULE_TYPE_BY_SOURCE_TYPE_SPECS:
            if getattr(detection.molecule, type_spec["attribute"]) is True:
                overall_type_counts[type_spec["key"]] += 1

        credited_source_keys = {
            _source_type_category(source.type)
            for source in detection.sources
        }
        for source_key in credited_source_keys:
            if source_key not in source_keys:
                continue
            source_counts[source_key] += 1
            for type_spec in MOLECULE_TYPE_BY_SOURCE_TYPE_SPECS:
                if getattr(detection.molecule, type_spec["attribute"]) is True:
                    counts_by_source[source_key][type_spec["key"]] += 1

    categories = []
    for source_spec in MOLECULE_TYPE_BY_SOURCE_SOURCE_SPECS:
        counts_by_type = {
            type_spec["key"]: counts_by_source[source_spec["key"]][type_spec["key"]]
            for type_spec in MOLECULE_TYPE_BY_SOURCE_TYPE_SPECS
        }
        categories.append(
            MoleculeTypesForSource(
                key=source_spec["key"],
                label=source_spec["label"],
                source_count=source_counts[source_spec["key"]],
                counts_by_type=counts_by_type,
            )
        )

    return MoleculeTypeBySourceData(
        categories=tuple(categories),
        molecule_count=len(first_detections),
        overall_type_counts={
            type_spec["key"]: overall_type_counts[type_spec["key"]]
            for type_spec in MOLECULE_TYPE_BY_SOURCE_TYPE_SPECS
        },
    )


def du_by_source_type_data(
    view: CensusView,
    detection_type: str = "ISM/CSM",
    *,
    include_tentative: bool = False,
    include_disputed: bool = False,
    include_isotopologues: bool = False,
    include_fullerenes: bool = False,
) -> DUBySourceTypeData:
    """Return DU values grouped by generalized first-detection source type.

    Each molecule is credited at most once per source category, matching the
    legacy convention used for source-type trend figures. Fullerenes are
    excluded by default because their extreme DU values dominate this
    comparison.
    """
    first_detections = first_detections_by_molecule(
        view,
        detection_type,
        include_tentative=include_tentative,
        include_disputed=include_disputed,
        include_isotopologues=include_isotopologues,
    )
    if not first_detections:
        raise ValueError("Cannot build a DU-by-source plot with no detections.")

    source_keys = {
        spec["key"]
        for spec in DU_BY_SOURCE_TYPE_SPECS
    }
    values_by_source = {
        source_key: []
        for source_key in source_keys
    }
    included_molecule_labels = set()

    for detection in first_detections.values():
        molecule = detection.molecule
        if molecule.fullerene and not include_fullerenes:
            continue
        if molecule.du is None:
            continue

        credited_source_keys = {
            _source_type_category(source.type)
            for source in detection.sources
        }
        used_source_keys = credited_source_keys & source_keys
        if not used_source_keys:
            continue

        included_molecule_labels.add(molecule.label)
        for source_key in used_source_keys:
            values_by_source[source_key].append(float(molecule.du))

    if not any(values_by_source.values()):
        raise ValueError("No DU values were available for the requested view.")

    categories = []
    for spec in DU_BY_SOURCE_TYPE_SPECS:
        categories.append(
            DUBySourceTypeCategory(
                key=spec["key"],
                label=spec["label"],
                color=spec["color"],
                values=tuple(values_by_source[spec["key"]]),
            )
        )

    return DUBySourceTypeData(
        categories=tuple(categories),
        molecule_count=len(included_molecule_labels),
    )


def relative_du_by_source_type_data(
    view: CensusView,
    detection_type: str = "ISM/CSM",
    *,
    include_tentative: bool = False,
    include_disputed: bool = False,
    include_isotopologues: bool = False,
    include_fullerenes: bool = False,
) -> RelativeDUBySourceTypeData:
    """Return relative-DU values grouped by first-detection source type.

    Relative DU is defined as ``du / maxdu`` using the same formula-domain
    convention as the legacy code. Each molecule is credited at most once per
    source category.
    """
    first_detections = first_detections_by_molecule(
        view,
        detection_type,
        include_tentative=include_tentative,
        include_disputed=include_disputed,
        include_isotopologues=include_isotopologues,
    )
    if not first_detections:
        raise ValueError(
            "Cannot build a relative-DU-by-source plot with no detections."
        )

    source_keys = {
        spec["key"]
        for spec in DU_BY_SOURCE_TYPE_SPECS
    }
    values_by_source = {
        source_key: []
        for source_key in source_keys
    }
    included_molecule_labels = set()

    for detection in first_detections.values():
        molecule = detection.molecule
        if molecule.fullerene and not include_fullerenes:
            continue
        if molecule.du is None or molecule.maxdu is None or molecule.maxdu == 0:
            continue

        credited_source_keys = {
            _source_type_category(source.type)
            for source in detection.sources
        }
        used_source_keys = credited_source_keys & source_keys
        if not used_source_keys:
            continue

        included_molecule_labels.add(molecule.label)
        value = float(molecule.du / molecule.maxdu)
        for source_key in used_source_keys:
            values_by_source[source_key].append(value)

    if not any(values_by_source.values()):
        raise ValueError(
            "No relative-DU values were available for the requested view."
        )

    categories = []
    for spec in DU_BY_SOURCE_TYPE_SPECS:
        categories.append(
            RelativeDUBySourceTypeCategory(
                key=spec["key"],
                label=spec["label"],
                color=spec["color"],
                values=tuple(values_by_source[spec["key"]]),
            )
        )

    return RelativeDUBySourceTypeData(
        categories=tuple(categories),
        molecule_count=len(included_molecule_labels),
    )


def mass_by_source_type_data(
    view: CensusView,
    detection_type: str = "ISM/CSM",
    *,
    include_tentative: bool = False,
    include_disputed: bool = False,
    include_isotopologues: bool = False,
    include_fullerenes: bool = False,
) -> MassBySourceTypeData:
    """Return molecular masses grouped by generalized first-detection source.

    Each molecule is credited at most once per source category, matching the
    legacy source-type figures. Fullerenes are excluded by default because
    their masses dominate this comparison.
    """
    first_detections = first_detections_by_molecule(
        view,
        detection_type,
        include_tentative=include_tentative,
        include_disputed=include_disputed,
        include_isotopologues=include_isotopologues,
    )
    if not first_detections:
        raise ValueError(
            "Cannot build a mass-by-source plot with no detections."
        )

    source_keys = {
        spec["key"]
        for spec in DU_BY_SOURCE_TYPE_SPECS
    }
    masses_by_source = {
        source_key: []
        for source_key in source_keys
    }
    included_molecule_labels = set()
    all_masses = []

    for detection in first_detections.values():
        molecule = detection.molecule
        if molecule.fullerene and not include_fullerenes:
            continue
        if molecule.mass is None:
            continue

        molecule_mass = float(molecule.mass)
        all_masses.append(molecule_mass)
        credited_source_keys = {
            _source_type_category(source.type)
            for source in detection.sources
        }
        used_source_keys = credited_source_keys & source_keys
        if not used_source_keys:
            continue

        included_molecule_labels.add(molecule.label)
        for source_key in used_source_keys:
            masses_by_source[source_key].append(molecule_mass)

    if not all_masses:
        raise ValueError("No molecular masses were available for the requested view.")
    if not any(masses_by_source.values()):
        raise ValueError(
            "No molecular masses had source categories for the requested view."
        )

    categories = []
    for spec in DU_BY_SOURCE_TYPE_SPECS:
        categories.append(
            MassBySourceTypeCategory(
                key=spec["key"],
                label=spec["label"],
                color=spec["color"],
                masses=tuple(masses_by_source[spec["key"]]),
            )
        )

    return MassBySourceTypeData(
        categories=tuple(categories),
        molecule_count=len(included_molecule_labels),
        mass_range=(min(all_masses), max(all_masses)),
    )


def wavelength_by_source_type_data(
    view: CensusView,
    detection_type: str = "ISM/CSM",
    *,
    include_tentative: bool = False,
    include_disputed: bool = False,
    include_isotopologues: bool = False,
    include_fullerenes: bool = True,
) -> WavelengthBySourceTypeData:
    """Return first-detection wavelength counts by generalized source type.

    Each molecule is credited at most once per source category, matching the
    legacy source/wavelength pie-grid convention. Within a credited source
    category, every wavelength listed for the first detection is counted.
    """
    first_detections = first_detections_by_molecule(
        view,
        detection_type,
        include_tentative=include_tentative,
        include_disputed=include_disputed,
        include_isotopologues=include_isotopologues,
    )
    if not first_detections:
        raise ValueError(
            "Cannot build wavelength-by-source data with no detections."
        )

    source_keys = {
        spec["key"]
        for spec in DU_BY_SOURCE_TYPE_SPECS
    }
    counts_by_source = {
        source_key: {
            wavelength: 0
            for wavelength in WAVES_BY_SOURCE_TYPE_WAVELENGTHS
        }
        for source_key in source_keys
    }
    included_molecule_labels = set()

    for detection in first_detections.values():
        molecule = detection.molecule
        if molecule.fullerene and not include_fullerenes:
            continue

        credited_source_keys = {
            _source_type_category(source.type)
            for source in detection.sources
        }
        used_source_keys = credited_source_keys & source_keys
        if not used_source_keys:
            continue

        credited_wavelengths = [
            wavelength
            for wavelength in WAVES_BY_SOURCE_TYPE_WAVELENGTHS
            if wavelength in detection.wavelengths
        ]
        if not credited_wavelengths:
            continue

        included_molecule_labels.add(molecule.label)
        for source_key in used_source_keys:
            for wavelength in credited_wavelengths:
                counts_by_source[source_key][wavelength] += 1

    if not any(sum(counts.values()) for counts in counts_by_source.values()):
        raise ValueError(
            "No wavelength/source counts were available for the requested view."
        )

    categories = []
    for spec in DU_BY_SOURCE_TYPE_SPECS:
        categories.append(
            WavelengthBySourceTypeCategory(
                key=spec["key"],
                label=spec["label"],
                counts_by_wavelength=counts_by_source[spec["key"]],
            )
        )

    return WavelengthBySourceTypeData(
        categories=tuple(categories),
        molecule_count=len(included_molecule_labels),
    )


def mass_by_wavelength_data(
    view: CensusView,
    detection_type: str = "ISM/CSM",
    *,
    include_tentative: bool = False,
    include_disputed: bool = False,
    include_isotopologues: bool = False,
    include_fullerenes: bool = True,
) -> MassByWavelengthData:
    """Return molecular masses grouped by first-detection wavelength."""
    first_detections = first_detections_by_molecule(
        view,
        detection_type,
        include_tentative=include_tentative,
        include_disputed=include_disputed,
        include_isotopologues=include_isotopologues,
    )
    if not first_detections:
        raise ValueError("Cannot build mass-by-wavelength data with no detections.")

    masses_by_wave = defaultdict(list)
    uv_vis_by_label = {}
    for detection in first_detections.values():
        molecule = detection.molecule
        if not include_fullerenes and molecule.fullerene:
            continue
        for wavelength in ("cm", "mm", "sub-mm", "IR", "UV", "Vis"):
            if wavelength in detection.wavelengths:
                masses_by_wave[wavelength].append(molecule.mass)
        if "UV" in detection.wavelengths or "Vis" in detection.wavelengths:
            uv_vis_by_label[molecule.label] = molecule.mass

    masses_by_wave["UV-Vis"] = list(uv_vis_by_label.values())

    series = []
    for spec in MASS_BY_WAVELENGTH_SPECS:
        key = spec["key"]
        series.append(
            MassByWavelengthSeries(
                key=key,
                label=spec["label"],
                color=spec["color"],
                masses=tuple(masses_by_wave[key]),
                annotation_xy=spec["annotation_xy"],
                truncate_at_data_max=spec.get("truncate_at_data_max", False),
                zorder=spec.get("zorder"),
            )
        )

    return MassByWavelengthData(
        series=tuple(series),
        molecule_count=len(first_detections),
    )


def molecules_by_wavelength_atoms_data(
    view: CensusView,
    detection_type: str = "ISM/CSM",
    *,
    include_tentative: bool = False,
    include_disputed: bool = False,
    include_isotopologues: bool = False,
    include_fullerenes: bool = False,
) -> MoleculesByWavelengthAtomsData:
    """Return molecule atom counts grouped by first-detection wavelength."""
    first_detections = first_detections_by_molecule(
        view,
        detection_type,
        include_tentative=include_tentative,
        include_disputed=include_disputed,
        include_isotopologues=include_isotopologues,
    )
    if not first_detections:
        raise ValueError("Cannot build wavelength atom-count data with no detections.")

    atom_counts_by_wave = defaultdict(list)
    included_molecules = set()
    for detection in first_detections.values():
        molecule = detection.molecule
        if not include_fullerenes and molecule.fullerene:
            continue
        for wavelength in ("UV", "Vis", "IR", "sub-mm", "mm", "cm"):
            if wavelength in detection.wavelengths:
                atom_counts_by_wave[wavelength].append(molecule.natoms)
                included_molecules.add(molecule.label)

    series = []
    for spec in MOLECULES_BY_WAVELENGTH_ATOMS_SPECS:
        key = spec["key"]
        series.append(
            MoleculesByWavelengthAtomsSeries(
                key=key,
                label=spec["label"],
                plot=spec["plot"],
                atom_counts=tuple(atom_counts_by_wave[key]),
            )
        )

    return MoleculesByWavelengthAtomsData(
        series=tuple(series),
        molecule_count=len(included_molecules),
    )


def detection_rate_by_atoms_data(
    view: CensusView,
    detection_type: str = "ISM/CSM",
    *,
    end_year: int | None = None,
    include_tentative: bool = False,
    include_disputed: bool = False,
    include_isotopologues: bool = False,
    categories: tuple[dict[str, object], ...] = DETECTION_RATE_BY_ATOMS_CATEGORIES,
) -> DetectionRateByAtomsData:
    """Return average detections/year by atom-count category.

    The rate is the number of first detections in a category divided by the
    number of years from that category's first detection through ``end_year``.
    This matches the legacy ``det_per_year_per_atom`` figure definition.
    """
    if end_year is None:
        end_year = view.census_year if not view.is_current else date.today().year

    first_detections = first_detections_by_molecule(
        view,
        detection_type,
        include_tentative=include_tentative,
        include_disputed=include_disputed,
        include_isotopologues=include_isotopologues,
    )
    if not first_detections:
        raise ValueError("Cannot build a rate-by-atoms plot with no detections.")

    points = []
    for category_spec in categories:
        category = category_spec["category"]
        category_years = [
            detection.year
            for detection in first_detections.values()
            if _detection_matches_atom_category(detection, category)
        ]
        if category_years:
            first_year = min(category_years)
            count = len(category_years)
            rate = count / (end_year - first_year + 1)
        else:
            first_year = None
            count = 0
            rate = np.nan

        points.append(
            DetectionRateByAtomsPoint(
                category=category,
                label=str(category_spec["label"]),
                x_position=float(category_spec["x"]),
                count=count,
                first_year=first_year,
                rate=float(rate),
            )
        )

    return DetectionRateByAtomsData(
        end_year=end_year,
        points=tuple(points),
    )


def facility_share_data(
    view: CensusView,
    detection_type: str = "ISM/CSM",
    *,
    end_year: int | None = None,
    top_n: int = 9,
    facility_nicks: tuple[str, ...] | None = None,
    include_tentative: bool = False,
    include_disputed: bool = False,
    include_isotopologues: bool = False,
    use_legacy_2021_selection: bool = True,
) -> FacilityShareData:
    """Return facility-share data for first molecule detections.

    Facilities are ranked by their number of first detections, then displayed
    by descending share of all first detections that occurred during their
    operational window. The 2021 published figure used a facility set that is
    preserved when ``use_legacy_2021_selection`` is true.
    """
    if end_year is None:
        end_year = view.census_year if not view.is_current else date.today().year

    first_detections = list(
        first_detections_by_molecule(
            view,
            detection_type,
            include_tentative=include_tentative,
            include_disputed=include_disputed,
            include_isotopologues=include_isotopologues,
        ).values()
    )
    if not first_detections:
        raise ValueError("Cannot build a facility-share plot with no detections.")

    detection_counts = defaultdict(int)
    for detection in first_detections:
        for telescope in detection.telescopes:
            detection_counts[telescope.nick] += 1

    if facility_nicks is None:
        if (
            use_legacy_2021_selection
            and not view.is_current
            and view.census_year == 2021
        ):
            facility_nicks = LEGACY_2021_FACILITY_SHARE_NICKS
        else:
            ranked_nicks = sorted(
                detection_counts,
                key=lambda nick: (
                    -detection_counts[nick],
                    view.db.telescopes[nick].shortname,
                ),
            )
            facility_nicks = tuple(ranked_nicks[:top_n])

    all_detection_years = [detection.year for detection in first_detections]
    facilities = []
    for nick in facility_nicks:
        telescope = view.db.telescopes[nick]
        if telescope.built is None:
            continue
        active_at_view = (
            telescope.decommissioned is None
            or telescope.decommissioned > end_year
        )
        facility_end_year = (
            telescope.decommissioned
            if telescope.decommissioned is not None
            and telescope.decommissioned <= end_year
            else end_year
        )
        total_window_detections = sum(
            telescope.built <= detection_year <= facility_end_year
            for detection_year in all_detection_years
        )
        facilities.append(
            FacilityShare(
                nick=nick,
                label=telescope.latex_name or telescope.shortname or telescope.name,
                start_year=telescope.built,
                end_year=facility_end_year,
                detection_count=detection_counts[nick],
                total_window_detections=total_window_detections,
                active_at_view=active_at_view,
            )
        )

    return FacilityShareData(
        end_year=end_year,
        facilities=tuple(
            sorted(
                facilities,
                key=lambda facility: (
                    facility.fraction,
                    facility.detection_count,
                    facility.label,
                ),
                reverse=True,
            )[:top_n]
        ),
    )


def scopes_by_year_data(
    view: CensusView,
    detection_type: str = "ISM/CSM",
    *,
    start_year: int = 1965,
    end_year: int | None = None,
    min_detections: int = 10,
    include_tentative: bool = False,
    include_disputed: bool = False,
    include_isotopologues: bool = False,
    cutoffs: dict[str, int] | None = None,
    colors_by_shortname: dict[str, str] | None = None,
) -> ScopesByYearData:
    """Return cumulative facility first-detection contributions by year."""
    if end_year is None:
        end_year = view.census_year if not view.is_current else date.today().year
    if end_year < start_year:
        raise ValueError("end_year must be greater than or equal to start_year.")

    if cutoffs is None:
        cutoffs = SCOPES_BY_YEAR_CUTOFFS
    if colors_by_shortname is None:
        colors_by_shortname = LEGACY_SCOPE_COLORS_BY_SHORTNAME

    first_detections = list(
        first_detections_by_molecule(
            view,
            detection_type,
            include_tentative=include_tentative,
            include_disputed=include_disputed,
            include_isotopologues=include_isotopologues,
        ).values()
    )
    if not first_detections:
        raise ValueError("Cannot build a scopes-by-year plot with no detections.")

    years = np.arange(start_year, end_year + 1)
    detections_by_telescope = defaultdict(list)
    for detection in first_detections:
        for telescope in detection.telescopes:
            detections_by_telescope[telescope.nick].append(detection)

    series = []
    for nick, detections in detections_by_telescope.items():
        if len(detections) < min_detections:
            continue

        telescope = view.db.telescopes[nick]
        if telescope.built is None:
            continue

        detection_years = [detection.year for detection in detections]
        counts = np.array(
            [sum(detection_year <= year for detection_year in detection_years) for year in years]
        )
        label = telescope.latex_name or telescope.shortname or telescope.name
        if telescope.shortname in cutoffs:
            fit_stop_year = min(cutoffs[telescope.shortname], end_year)
            fit_mask = (years >= telescope.built) & (years < fit_stop_year)
        else:
            fit_stop_year = end_year
            fit_mask = (years >= telescope.built) & (years <= fit_stop_year)
        if fit_mask.sum() < 2:
            rate = 0.0
        else:
            rate = (
                np.polynomial.polynomial.Polynomial.fit(
                    years[fit_mask],
                    counts[fit_mask],
                    1,
                )
                .convert()
                .coef[1]
            )

        series.append(
            ScopeDetectionSeries(
                nick=nick,
                label=label,
                start_year=telescope.built,
                fit_stop_year=fit_stop_year,
                color=colors_by_shortname.get(telescope.shortname, ASTROMOL_BLUE),
                years=years,
                counts=counts,
                detection_count=len(detections),
                rate=float(rate),
            )
        )

    return ScopesByYearData(
        years=years,
        series=tuple(
            sorted(
                series,
                key=lambda item: (-item.detection_count, item.label),
            )
        ),
        end_year=end_year,
    )


def cumulative_counts_by_year(
    first_detection_years: dict[str, int],
    *,
    start_year: int | None = None,
    end_year: int | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Return cumulative detection counts for each year in a range."""
    if not first_detection_years:
        raise ValueError("Cannot build a cumulative series with no detections.")

    detection_years = list(first_detection_years.values())
    if start_year is None:
        start_year = min(detection_years)
    if end_year is None:
        end_year = max(detection_years)
    if end_year < start_year:
        raise ValueError("end_year must be greater than or equal to start_year.")

    years = np.arange(start_year, end_year + 1)
    counts = np.array(
        [
            sum(detection_year <= year for detection_year in detection_years)
            for year in years
        ]
    )
    return years, counts


def linear_cumulative_rate(
    years: np.ndarray,
    counts: np.ndarray,
    *,
    start_year: int,
    stop_year: int | None = None,
) -> float:
    """Fit a linear cumulative-detection rate over a year range.

    ``stop_year`` is included in the fit. This makes plot annotations read as
    closed date ranges, such as 1968-2005 and 2005-2021.
    """
    if start_year not in years:
        raise ValueError(f"start_year {start_year} is outside the data range.")
    if stop_year is not None and stop_year not in years:
        raise ValueError(f"stop_year {stop_year} is outside the data range.")

    start_index = int(np.argwhere(years == start_year)[0][0])
    stop_index = (
        int(np.argwhere(years == stop_year)[0][0]) + 1
        if stop_year is not None
        else None
    )
    fit_years = years[start_index:stop_index]
    fit_counts = counts[start_index:stop_index]
    if len(fit_years) < 2:
        raise ValueError("At least two years are required for a linear fit.")

    return float(
        np.polynomial.polynomial.Polynomial.fit(
            fit_years,
            fit_counts,
            1,
        ).convert().coef[1]
    )


def cumulative_detection_trend_ranges(
    view: CensusView,
    end_year: int,
) -> tuple[tuple[str, int, int | None], ...]:
    """Return census-aware trend ranges for cumulative-detection plots."""
    if view.is_current:
        return (
            ("1968-2005", 1968, 2005),
            ("2005-2021", 2005, 2021),
            ("2021-Present", 2021, None),
        )

    census_year = view.census_year
    if census_year is None:
        raise ValueError("Census trend ranges require a census year.")

    if census_year <= 2021:
        return (
            ("1968-2005", 1968, 2005),
            (f"2005-{end_year}", 2005, end_year),
        )

    return (
        ("1968-2005", 1968, 2005),
        ("2005-2021", 2005, 2021),
        (f"2021-{end_year}", 2021, end_year),
    )


def cumulative_detection_data(
    view: CensusView,
    detection_type: str = "ISM/CSM",
    *,
    start_year: int | None = None,
    end_year: int | None = None,
    include_tentative: bool = False,
    include_disputed: bool = False,
    include_isotopologues: bool = False,
) -> CumulativeDetectionData:
    """Return cumulative-detection data for a census view.

    By default, the data represent secure, non-isotopologue ISM/CSM molecule
    detections. For census views, the default end year is the census boundary;
    for current views, it is the current calendar year.
    """
    if end_year is None:
        end_year = view.census_year if not view.is_current else date.today().year

    first_years = first_detection_years_by_molecule(
        view,
        detection_type,
        include_tentative=include_tentative,
        include_disputed=include_disputed,
        include_isotopologues=include_isotopologues,
    )
    years, counts = cumulative_counts_by_year(
        first_years,
        start_year=start_year,
        end_year=end_year,
    )

    trends = []
    for label, trend_start, trend_stop in cumulative_detection_trend_ranges(
        view,
        end_year,
    ):
        if trend_start < years[0] or trend_start > years[-1]:
            continue
        if trend_stop is not None and trend_stop > years[-1]:
            continue
        trends.append(
            DetectionTrend(
                label=label,
                start_year=trend_start,
                stop_year=trend_stop,
                slope=linear_cumulative_rate(
                    years,
                    counts,
                    start_year=trend_start,
                    stop_year=trend_stop,
                ),
            )
        )

    return CumulativeDetectionData(
        years=years,
        counts=counts,
        first_detection_years=first_years,
        trends=tuple(trends),
        show_facility_markers=(
            not view.is_current and view.census_year in {2018, 2021}
        ),
    )


def _pyplot():
    """Import matplotlib only when a plotting function is used."""
    import matplotlib.pyplot as plt

    return plt


def _style_manuscript_axes(
    ax,
    *,
    xlabel: str,
    ylabel: str,
    text_size: int = FIGURE_TEXT_SIZE,
) -> None:
    """Apply the common manuscript style used by astromol figures."""
    ax.set_xlabel(xlabel, fontsize=text_size)
    ax.set_ylabel(ylabel, fontsize=text_size)
    ax.tick_params(
        axis="x",
        which="both",
        direction="in",
        length=FIGURE_TICK_LENGTH,
        width=FIGURE_TICK_WIDTH,
        labelsize=text_size,
    )
    ax.tick_params(
        axis="y",
        which="both",
        direction="in",
        length=FIGURE_TICK_LENGTH,
        width=FIGURE_TICK_WIDTH,
        labelsize=text_size,
    )
    ax.yaxis.set_ticks_position("both")
    ax.xaxis.set_ticks_position("both")


def _finalize_manuscript_figure(figure, ax) -> None:
    """Use a fixed axes rectangle so paired figures align visually."""
    ax.set_position(FIGURE_AXES_BOUNDS)


def _padded_count_axis_limit(max_count: int, padding_fraction: float = 0.12) -> int:
    """Return a rounded y-axis limit with room above a cumulative curve."""
    if max_count < 1:
        return 1

    padded_count = max_count * (1 + padding_fraction)
    return int(np.ceil(padded_count / 10) * 10)


def plot_cumulative_detections(
    data: CumulativeDetectionData,
    *,
    ax=None,
    color: str = ASTROMOL_BLUE,
    linewidth: float = 4,
    annotate_trends: bool = True,
    annotate_total: bool = True,
    annotate_facilities: bool | None = None,
    text_size: int = 24,
):
    """Plot cumulative detections and return ``(figure, axes)``."""
    plt = _pyplot()
    if annotate_facilities is None:
        annotate_facilities = data.show_facility_markers

    if ax is None:
        figure, ax = plt.subplots(
            num="Cumulative Detections",
            figsize=FIGURE_SIZE,
        )
    else:
        figure = ax.figure

    ax.plot(data.years, data.counts, color=color, linewidth=linewidth)
    ax.set_ylim(0, _padded_count_axis_limit(data.total))
    _style_manuscript_axes(
        ax,
        xlabel="Year",
        ylabel="Cumulative Number of Detected Molecules",
        text_size=text_size,
    )

    if annotate_trends:
        trend_lines = [
            f"{trend.label}: {trend.slope:.1f} detections/year"
            for trend in data.trends
        ]
        ax.annotate(
            "\n".join(trend_lines),
            xy=(0.05, 0.95),
            xycoords="axes fraction",
            ha="left",
            va="top",
            size=text_size,
        )

    if annotate_total:
        ax.annotate(
            f"Total: {data.total}",
            xy=(0.95, 0.95),
            xycoords="axes fraction",
            va="top",
            ha="right",
            size=text_size,
        )

    if annotate_facilities:
        arrowprops = {
            "arrowstyle": "-|>",
            "connectionstyle": "arc3",
            "facecolor": "black",
        }
        counts_by_year = dict(zip(data.years, data.counts))
        for marker in CUMULATIVE_DETECTION_FACILITY_MARKERS:
            marker_year = marker["year"]
            if marker_year not in counts_by_year:
                continue
            ax.annotate(
                marker["label"],
                xy=(
                    marker_year,
                    counts_by_year[marker_year] + marker["count_offset"],
                ),
                xycoords="data",
                xytext=marker["text_offset"],
                textcoords="offset points",
                rotation=marker["rotation"],
                arrowprops=arrowprops,
                va=marker["vertical_alignment"],
                ha=marker["horizontal_alignment"],
                size=text_size,
            )

    _finalize_manuscript_figure(figure, ax)
    return figure, ax


def write_cumulative_detections_plot(
    data: CumulativeDetectionData,
    output_path: str | Path,
    *,
    file_format: str | None = None,
) -> Path:
    """Write the cumulative-detections plot and return the output path."""
    output_path = Path(output_path)
    if file_format is None:
        file_format = output_path.suffix.lstrip(".") or "pdf"

    figure, _ = plot_cumulative_detections(data)
    figure.savefig(
        output_path,
        format=file_format,
        transparent=True,
    )
    _pyplot().close(figure)
    return output_path


def _annotate_atom_count_legend(
    ax,
    series: tuple[CumulativeByAtomsSeries, ...],
    *,
    legend_text_size: int = FIGURE_LEGEND_TEXT_SIZE,
    legend_x: float = 0.10,
    legend_y: float = 0.92,
    legend_step: float = 0.04,
) -> None:
    """Annotate atom-count category labels inside a plot, legacy style."""
    for index, series_item in enumerate(series):
        ax.annotate(
            series_item.label,
            xy=(legend_x, legend_y - (index * legend_step)),
            xycoords="axes fraction",
            color=series_item.color,
            size=legend_text_size,
            weight="bold",
            va="top",
            ha="left",
        )


def plot_cumulative_by_atoms(
    data: CumulativeByAtomsData,
    *,
    ax=None,
    linewidth: float = 2.5,
    text_size: int = FIGURE_TEXT_SIZE,
    legend_text_size: int = FIGURE_LEGEND_TEXT_SIZE,
    legend_x: float = 0.10,
    legend_y: float = 0.92,
    legend_step: float = 0.04,
    y_axis_padding_fraction: float = 0.05,
):
    """Plot cumulative ISM/CSM detections grouped by atom-count category."""
    plt = _pyplot()
    if ax is None:
        figure, ax = plt.subplots(
            num="Cumulative Detections By Atoms",
            figsize=FIGURE_SIZE,
        )
    else:
        figure = ax.figure

    for series in data.series:
        ax.plot(
            series.years,
            series.counts,
            color=series.color,
            linewidth=linewidth,
        )

    _style_manuscript_axes(
        ax,
        xlabel="Year",
        ylabel="Cumulative Number of Detected Molecules",
        text_size=text_size,
    )
    ax.set_xlim(data.start_year, data.end_year + 3)
    ax.set_ylim(-1, data.max_count * (1 + y_axis_padding_fraction))
    ax.set_xticks(
        [
            tick
            for tick in range(1940, data.end_year + 20, 20)
            if data.start_year <= tick <= data.end_year + 3
        ]
    )

    _annotate_atom_count_legend(
        ax,
        data.series,
        legend_text_size=legend_text_size,
        legend_x=legend_x,
        legend_y=legend_y,
        legend_step=legend_step,
    )

    _finalize_manuscript_figure(figure, ax)
    return figure, ax


def write_cumulative_by_atoms_plot(
    data: CumulativeByAtomsData,
    output_path: str | Path,
    *,
    file_format: str | None = None,
) -> Path:
    """Write the cumulative-by-atoms plot and return the output path."""
    output_path = Path(output_path)
    if file_format is None:
        file_format = output_path.suffix.lstrip(".") or "pdf"

    figure, _ = plot_cumulative_by_atoms(data)
    figure.savefig(
        output_path,
        format=file_format,
        transparent=True,
    )
    _pyplot().close(figure)
    return output_path


def plot_stacked_cumulative_by_atoms(
    data: CumulativeByAtomsData,
    *,
    ax=None,
    text_size: int = FIGURE_TEXT_SIZE,
    legend_text_size: int = FIGURE_LEGEND_TEXT_SIZE,
    legend_x: float = 0.10,
    legend_y: float = 0.92,
    legend_step: float = 0.04,
    stack_alpha: float = 0.88,
    y_axis_padding_fraction: float = 0.04,
):
    """Plot cumulative ISM/CSM detections as a stacked atom-count history."""
    plt = _pyplot()
    if ax is None:
        figure, ax = plt.subplots(
            num="Stacked Cumulative Detections By Atoms",
            figsize=FIGURE_SIZE,
        )
    else:
        figure = ax.figure

    active_series = [
        series for series in data.series if series.final_count > 0
    ]
    ax.stackplot(
        data.years,
        [series.counts for series in active_series],
        colors=[series.color for series in active_series],
        alpha=stack_alpha,
        linewidth=0.3,
        edgecolor="white",
    )

    _style_manuscript_axes(
        ax,
        xlabel="Year",
        ylabel="Cumulative Number of Detected Molecules",
        text_size=text_size,
    )
    ax.set_xlim(data.start_year, data.end_year + 3)
    ax.set_ylim(0, data.total * (1 + y_axis_padding_fraction))
    ax.set_xticks(
        [
            tick
            for tick in range(1940, data.end_year + 20, 20)
            if data.start_year <= tick <= data.end_year + 3
        ]
    )

    _annotate_atom_count_legend(
        ax,
        data.series,
        legend_text_size=legend_text_size,
        legend_x=legend_x,
        legend_y=legend_y,
        legend_step=legend_step,
    )

    _finalize_manuscript_figure(figure, ax)
    return figure, ax


def write_stacked_cumulative_by_atoms_plot(
    data: CumulativeByAtomsData,
    output_path: str | Path,
    *,
    file_format: str | None = None,
) -> Path:
    """Write the stacked cumulative-by-atoms plot and return the output path."""
    output_path = Path(output_path)
    if file_format is None:
        file_format = output_path.suffix.lstrip(".") or "pdf"

    figure, _ = plot_stacked_cumulative_by_atoms(data)
    figure.savefig(
        output_path,
        format=file_format,
        transparent=True,
    )
    _pyplot().close(figure)
    return output_path


def _periodic_heatmap_color(count: int, max_count: int) -> str:
    """Return the legacy yellow-to-red heatmap color for one element count."""
    from matplotlib.colors import to_hex, to_rgb

    if count <= 0:
        return "white"
    if max_count <= 1:
        fraction = 1.0
    else:
        fraction = (count - 1) / (max_count - 1)
    start = np.array(to_rgb(PERIODIC_HEATMAP_START_COLOR))
    stop = np.array(to_rgb(PERIODIC_HEATMAP_STOP_COLOR))
    return to_hex(start + (stop - start) * fraction)


def _periodic_heatmap_colormap():
    """Return the yellow-to-red colormap used for astromol heatmaps."""
    from matplotlib.colors import LinearSegmentedColormap

    return LinearSegmentedColormap.from_list(
        "astromol_yellow_red",
        [PERIODIC_HEATMAP_START_COLOR, PERIODIC_HEATMAP_STOP_COLOR],
    )


def plot_periodic_heatmap(
    data: PeriodicHeatmapData,
    *,
    ax=None,
    figure_size: tuple[float, float] = (20, 9.5),
):
    """Plot a periodic-table heatmap of element occurrence in molecules."""
    plt = _pyplot()
    import matplotlib.patches as patches

    if ax is None:
        figure = plt.figure(num="Periodic Heatmap", figsize=figure_size)
        ax = figure.add_axes([0, 0, 1, 1])
    else:
        figure = ax.figure

    ax.set_xlim(0, 18)
    ax.set_ylim(0, 8)

    for cell in data.cells:
        x_position = cell.x_position
        y_position = cell.y_position
        rect = patches.Rectangle(
            (x_position, y_position),
            0.9,
            1.05,
            linewidth=1,
            edgecolor="black",
            facecolor=_periodic_heatmap_color(cell.count, data.max_count),
            alpha=0.5,
        )
        ax.add_patch(rect)

        if cell.count:
            ax.annotate(
                str(cell.count),
                xy=(x_position + 0.8, y_position + 0.95),
                xycoords="data",
                size=14,
                color="black",
                ha="right",
                va="top",
                weight="bold",
            )

        ax.annotate(
            str(cell.atomic_number),
            xy=(x_position + 0.1, y_position + 0.95),
            xycoords="data",
            size=14,
            color="black",
            ha="left",
            va="top",
        )
        ax.annotate(
            cell.symbol,
            xy=(x_position + 0.1, y_position + 0.70),
            xycoords="data",
            size=20,
            color="black",
            ha="left",
            va="top",
            weight="bold",
        )
        ax.annotate(
            f"{cell.atomic_mass:.3f}",
            xy=(x_position + 0.1, y_position + 0.42),
            xycoords="data",
            size=8,
            color="black",
            ha="left",
            va="top",
        )
        ax.annotate(
            cell.name,
            xy=(x_position + 0.1, y_position + 0.29),
            xycoords="data",
            size=8,
            color="black",
            ha="left",
            va="top",
        )

    ax.set_aspect("equal")
    ax.axis("off")
    return figure, ax


def write_periodic_heatmap(
    data: PeriodicHeatmapData,
    output_path: str | Path,
    *,
    file_format: str | None = None,
) -> Path:
    """Write the periodic-table heatmap and return the output path."""
    output_path = Path(output_path)
    if file_format is None:
        file_format = output_path.suffix.lstrip(".") or "pdf"

    figure, _ = plot_periodic_heatmap(data)
    figure.savefig(
        output_path,
        format=file_format,
        transparent=True,
        bbox_inches="tight",
        pad_inches=0,
    )
    _pyplot().close(figure)
    return output_path


def plot_du_histogram(
    data: DUHistogramData,
    *,
    ax=None,
    bins: tuple[float, ...] | None = None,
    text_size: int = FIGURE_TEXT_SIZE,
    annotation_size: int = FIGURE_LEGEND_TEXT_SIZE,
    facecolor: str = ASTROMOL_BLUE,
    edgecolor: str = "royalblue",
    alpha: float = 0.25,
):
    """Plot the legacy degree-of-unsaturation histogram."""
    plt = _pyplot()
    if bins is None:
        bins = _du_histogram_bins(data)

    if ax is None:
        figure, ax = plt.subplots(
            num="Degree of Unsaturation Histogram",
            figsize=FIGURE_SIZE,
        )
    else:
        figure = ax.figure

    bin_array = np.array(bins)
    counts, _, _ = ax.hist(
        data.values,
        bins=bin_array,
        facecolor=facecolor,
        alpha=alpha,
    )
    ax.hist(
        data.values,
        bins=bin_array,
        edgecolor=edgecolor,
        linewidth=1.5,
        fill=False,
    )

    ax.set_xlabel("Degree of Unsaturation", fontsize=text_size)
    ax.set_ylabel("# of Detected Molecules", fontsize=text_size)
    y_limit = max(35, int(np.ceil((max(counts) + 2) / 5) * 5))
    ax.set_ylim(0, y_limit)
    ax.tick_params(
        axis="x",
        which="both",
        direction="in",
        length=5,
        width=FIGURE_TICK_WIDTH,
        labelsize=text_size,
        top=True,
    )
    ax.tick_params(
        axis="y",
        which="both",
        direction="in",
        length=5,
        width=FIGURE_TICK_WIDTH,
        labelsize=text_size,
        right=True,
    )
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
    ax.set_title("")

    first_bin_count = counts[0] if len(counts) else 0
    ax.annotate(
        "CH$_4$, CH$_3$OH, CH$_3$Cl, ...",
        xy=(0, first_bin_count + 1),
        xycoords="data",
        rotation=90,
        size=annotation_size,
        ha="center",
        va="bottom",
    )

    hc11n_bin_index = 24
    hc11n_bin_count = counts[hc11n_bin_index] if len(counts) > hc11n_bin_index else 0
    ax.annotate(
        "HC$_{11}$N",
        xy=(12, hc11n_bin_count + 1),
        xycoords="data",
        rotation=90,
        size=annotation_size,
        ha="center",
        va="bottom",
    )

    if data.max_du > 12:
        max_du_formulas = [
            formula
            for formula, value in zip(data.formula_labels, data.values)
            if value == data.max_du
        ]
        max_du_bin_index = int(np.searchsorted(bin_array, data.max_du, side="right") - 1)
        max_du_bin_count = (
            counts[max_du_bin_index]
            if 0 <= max_du_bin_index < len(counts)
            else 0
        )
        ax.annotate(
            _matplotlib_formula_list(max_du_formulas),
            xy=(data.max_du, max_du_bin_count + 1),
            xycoords="data",
            rotation=90,
            size=annotation_size,
            ha="center",
            va="bottom",
        )

    figure.tight_layout()
    return figure, ax


def write_du_histogram(
    data: DUHistogramData,
    output_path: str | Path,
    *,
    file_format: str | None = None,
) -> Path:
    """Write the degree-of-unsaturation histogram and return the output path."""
    output_path = Path(output_path)
    if file_format is None:
        file_format = output_path.suffix.lstrip(".") or "pdf"

    figure, _ = plot_du_histogram(data)
    figure.savefig(
        output_path,
        format=file_format,
        transparent=True,
        bbox_inches="tight",
    )
    _pyplot().close(figure)
    return output_path


def plot_kappa_histogram(
    data: KappaHistogramData,
    *,
    ax=None,
    bins: int = 100,
    figure_size: tuple[float, float] = (8.0, 4.8),
    text_size: int = FIGURE_TEXT_SIZE,
    facecolor: str = ASTROMOL_BLUE,
    edgecolor: str = ASTROMOL_BLUE,
    alpha: float = 0.25,
    show_shape_guide: bool = True,
):
    """Plot the Ray asymmetry-parameter histogram."""
    plt = _pyplot()

    if ax is None:
        figure = plt.figure(num="Kappas", figsize=figure_size)
        ax = figure.add_subplot(111)
    else:
        figure = ax.figure

    ax.hist(data.values, bins, facecolor=facecolor, alpha=alpha)
    ax.hist(
        data.values,
        bins,
        edgecolor=edgecolor,
        linewidth=1.5,
        fill=False,
    )

    from matplotlib.ticker import ScalarFormatter

    ax.set_xlabel(r"$\kappa$", fontsize=text_size)
    ax.set_ylabel("# Molecules", fontsize=text_size)
    ax.set_xticks([-1.0, -0.5, 0.0, 0.5, 1.0])
    ax.set_yscale("log")
    ax.yaxis.set_major_formatter(ScalarFormatter())
    ax.tick_params(
        axis="x",
        which="both",
        direction="in",
        length=5,
        width=FIGURE_TICK_WIDTH,
        labelsize=text_size,
        colors="black",
    )
    ax.tick_params(
        axis="y",
        which="both",
        direction="in",
        length=5,
        width=FIGURE_TICK_WIDTH,
        labelsize=text_size,
        colors="black",
    )
    ax.xaxis.set_ticks_position("both")
    ax.yaxis.set_ticks_position("both")
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(FIGURE_TICK_WIDTH)
    if show_shape_guide:
        guide_text_size = text_size - 7
        y_guide = 1.06
        guide_transform = ax.get_xaxis_transform()
        ax.annotate(
            "",
            xy=(1.0, y_guide),
            xytext=(-1.0, y_guide),
            xycoords=guide_transform,
            textcoords=guide_transform,
            arrowprops={
                "arrowstyle": "<->",
                "color": "black",
                "linewidth": 1.0,
                "shrinkA": 0,
                "shrinkB": 0,
            },
            annotation_clip=False,
        )
        ax.text(
            -1.0,
            1.095,
            "prolate",
            transform=guide_transform,
            ha="left",
            va="bottom",
            fontsize=guide_text_size,
            clip_on=False,
        )
        ax.text(
            0.0,
            1.095,
            "asymmetric",
            transform=guide_transform,
            ha="center",
            va="bottom",
            fontsize=guide_text_size,
            clip_on=False,
        )
        ax.text(
            1.0,
            1.095,
            "oblate",
            transform=guide_transform,
            ha="right",
            va="bottom",
            fontsize=guide_text_size,
            clip_on=False,
        )
    ax.set_title("")
    return figure, ax


def write_kappa_histogram(
    data: KappaHistogramData,
    output_path: str | Path,
    *,
    bins: int = 100,
    file_format: str | None = None,
) -> Path:
    """Write the Ray asymmetry-parameter histogram and return the output path."""
    output_path = Path(output_path)
    if file_format is None:
        file_format = output_path.suffix.lstrip(".") or "pdf"

    figure, _ = plot_kappa_histogram(data, bins=bins)
    figure.savefig(
        output_path,
        format=file_format,
        transparent=True,
        bbox_inches="tight",
        pad_inches=0,
    )
    _pyplot().close(figure)
    return output_path


def _annotate_du_notable_formula(
    ax,
    *,
    data: DUHistogramData,
    value: float,
    count: int,
    annotation_size: int,
) -> None:
    """Annotate notable DU values using compact formula labels."""
    formulas = [
        formula
        for formula, du_value in zip(data.formula_labels, data.values)
        if du_value == value
    ]
    if not formulas:
        return
    ax.annotate(
        _matplotlib_formula_list(formulas),
        xy=(value, count + 1),
        xycoords="data",
        rotation=90,
        size=annotation_size,
        ha="center",
        va="bottom",
    )


def plot_du_bar_chart(
    data: DUHistogramData,
    *,
    ax=None,
    text_size: int = FIGURE_TEXT_SIZE,
    annotation_size: int = FIGURE_LEGEND_TEXT_SIZE,
    facecolor: str = ASTROMOL_BLUE,
    edgecolor: str = "black",
    alpha: float = 0.65,
    bar_width: float = 0.32,
    include_negative_du: bool = False,
):
    """Plot exact degree-of-unsaturation counts as a discrete bar chart."""
    plt = _pyplot()
    if ax is None:
        figure, ax = plt.subplots(
            num="Degree of Unsaturation Exact Counts",
            figsize=FIGURE_SIZE,
        )
    else:
        figure = ax.figure

    counts_by_value = data.value_counts(include_negative_du=include_negative_du)
    x_values = np.array(sorted(counts_by_value))
    y_values = np.array([counts_by_value[value] for value in x_values])

    ax.bar(
        x_values,
        y_values,
        width=bar_width,
        color=facecolor,
        alpha=alpha,
        edgecolor=edgecolor,
        linewidth=0.8,
    )

    ax.set_xlabel("Degree of Unsaturation", fontsize=text_size)
    ax.set_ylabel("# of Detected Molecules", fontsize=text_size)
    y_limit = int(np.ceil((max(y_values) + 6) / 5) * 5)
    ax.set_ylim(0, y_limit)
    ax.set_xlim(min(-0.9, min(x_values) - 0.4), max(14.7, max(x_values) + 0.4))
    ax.set_xticks(np.arange(0, max(15, int(np.ceil(max(x_values))) + 1), 2))
    ax.tick_params(
        axis="x",
        which="both",
        direction="in",
        length=5,
        width=FIGURE_TICK_WIDTH,
        labelsize=text_size,
        top=True,
    )
    ax.tick_params(
        axis="y",
        which="both",
        direction="in",
        length=5,
        width=FIGURE_TICK_WIDTH,
        labelsize=text_size,
        right=True,
    )
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
    ax.set_title("")

    if 0.0 in counts_by_value:
        ax.annotate(
            "CH$_4$, CH$_3$OH, CH$_3$Cl, ...",
            xy=(0, counts_by_value[0.0] + 1),
            xycoords="data",
            rotation=90,
            size=annotation_size,
            ha="center",
            va="bottom",
        )
    if 12.0 in counts_by_value:
        _annotate_du_notable_formula(
            ax,
            data=data,
            value=12.0,
            count=counts_by_value[12.0],
            annotation_size=annotation_size,
        )
    if data.max_du > 12:
        _annotate_du_notable_formula(
            ax,
            data=data,
            value=data.max_du,
            count=counts_by_value[data.max_du],
            annotation_size=annotation_size,
        )

    figure.tight_layout()
    return figure, ax


def write_du_bar_chart(
    data: DUHistogramData,
    output_path: str | Path,
    *,
    file_format: str | None = None,
) -> Path:
    """Write the exact-value degree-of-unsaturation bar chart."""
    output_path = Path(output_path)
    if file_format is None:
        file_format = output_path.suffix.lstrip(".") or "pdf"

    figure, _ = plot_du_bar_chart(data)
    figure.savefig(
        output_path,
        format=file_format,
        transparent=True,
        bbox_inches="tight",
    )
    _pyplot().close(figure)
    return output_path


def _type_pie_start_angle(fraction: float) -> float:
    """Return the legacy start angle for a molecule-type ring."""
    return -(90 - (360 - 360 * fraction) / 2)


def _ring_categories_by_order(categories, order: str):
    """Return ring categories in legacy or descending-count order."""
    if order == "legacy":
        return tuple(categories)
    if order == "count":
        return tuple(
            category
            for _, category in sorted(
                enumerate(categories),
                key=lambda indexed_category: (
                    -indexed_category[1].count,
                    indexed_category[0],
                ),
            )
        )
    raise ValueError("order must be 'legacy' or 'count'.")


def plot_type_pie_chart(
    data: MoleculeTypeData,
    *,
    ax=None,
    figure_size: tuple[float, float] = FIGURE_SIZE,
    ring_width: float = 0.10,
    ring_gap: float = 0.02,
    label_text_size: int = 13,
    label_y_offsets: tuple[float, ...] | None = None,
    percent_text_size: int = 11,
    remainder_color: str = "#EEEEEE",
    order: str = "count",
):
    """Plot the molecule-type concentric ring chart."""
    plt = _pyplot()

    if ax is None:
        figure = plt.figure(num="Type Pie Chart", figsize=figure_size)
        ax = figure.add_subplot(111)
    else:
        figure = ax.figure

    categories = _ring_categories_by_order(data.categories, order)

    for index, category in enumerate(categories):
        radius = 1 - index * (ring_width + ring_gap)
        ax.pie(
            [category.fraction, 1.0 - category.fraction],
            colors=[category.color, remainder_color],
            radius=radius,
            startangle=_type_pie_start_angle(category.fraction),
            wedgeprops={
                "width": ring_width,
                "edgecolor": "white",
                "linewidth": 1,
            },
        )

    label_y_positions = (0.110, 0.160, 0.205, 0.255, 0.305, 0.3575, 0.400)
    if label_y_offsets is None:
        label_y_offsets = (0.0,) * len(categories)
    if len(label_y_offsets) != len(categories):
        raise ValueError("label_y_offsets must match the number of type categories.")
    for category, y_position, y_offset in zip(
        categories,
        label_y_positions,
        label_y_offsets,
    ):
        ax.annotate(
            category.plural_label,
            xy=(0.5, y_position + y_offset),
            xycoords="axes fraction",
            color=category.color,
            ha="center",
            size=label_text_size,
            weight="bold",
        )

    for index, category in enumerate(categories):
        radius = 1 - index * (ring_width + ring_gap)
        radial_center = radius - ring_width / 2
        ax.annotate(
            f"{category.percent:.1f}%",
            xy=(radial_center, 0),
            xycoords="data",
            color=category.percent_color,
            ha="center",
            va="center",
            size=percent_text_size,
            rotation=-90,
            weight="bold",
        )

    ax.set_aspect("equal")
    ax.axis("off")
    figure.tight_layout()
    return figure, ax


def write_type_pie_chart(
    data: MoleculeTypeData,
    output_path: str | Path,
    *,
    file_format: str | None = None,
    order: str = "count",
) -> Path:
    """Write the molecule-type ring chart and return the output path."""
    output_path = Path(output_path)
    if file_format is None:
        file_format = output_path.suffix.lstrip(".") or "pdf"

    figure, _ = plot_type_pie_chart(data, order=order)
    figure.savefig(
        output_path,
        format=file_format,
        transparent=True,
        bbox_inches="tight",
    )
    _pyplot().close(figure)
    return output_path


def plot_source_pie_chart(
    data: SourceTypeData,
    *,
    ax=None,
    figure_size: tuple[float, float] = FIGURE_SIZE,
    ring_width: float = 0.10,
    ring_gap: float = 0.02,
    label_text_size: int = 14,
    percent_text_size: int = 12,
    remainder_color: str = "#EEEEEE",
    diffuse_cloud_label: str = "Diffuse Cloud",
    order: str = "count",
):
    """Plot the generalized source-type concentric ring chart."""
    plt = _pyplot()

    if ax is None:
        figure = plt.figure(num="Source Pie Chart", figsize=figure_size)
        ax = figure.add_subplot(111)
    else:
        figure = ax.figure

    categories = _ring_categories_by_order(data.categories, order)

    for index, category in enumerate(categories):
        radius = 1 - index * (ring_width + ring_gap)
        ax.pie(
            [category.fraction, 1.0 - category.fraction],
            colors=[category.color, remainder_color],
            radius=radius,
            startangle=_type_pie_start_angle(category.fraction),
            wedgeprops={
                "width": ring_width,
                "edgecolor": "white",
                "linewidth": 1,
            },
        )

    label_y_positions = (0.110, 0.160, 0.205, 0.255, 0.305)
    percent_y_positions = (0.870, 0.825, 0.775, 0.725, 0.680)
    for category, label_y, percent_y in zip(
        categories,
        label_y_positions,
        percent_y_positions,
    ):
        label = (
            diffuse_cloud_label
            if category.key == "diffuse_cloud"
            else category.label
        )
        ax.annotate(
            label,
            xy=(0.5, label_y),
            xycoords="axes fraction",
            color=category.color,
            ha="center",
            size=label_text_size,
            weight="bold",
        )
        ax.annotate(
            f"{category.percent:.1f}%",
            xy=(0.5, percent_y),
            xycoords="axes fraction",
            color="white",
            ha="center",
            size=percent_text_size,
            weight="bold",
        )

    ax.set_aspect("equal")
    ax.axis("off")
    figure.tight_layout()
    return figure, ax


def write_source_pie_chart(
    data: SourceTypeData,
    output_path: str | Path,
    *,
    file_format: str | None = None,
    diffuse_cloud_label: str = "Diffuse Cloud",
    pad_inches: float = -0.65,
    order: str = "count",
) -> Path:
    """Write the generalized source-type ring chart and return the path."""
    output_path = Path(output_path)
    if file_format is None:
        file_format = output_path.suffix.lstrip(".") or "pdf"

    figure, _ = plot_source_pie_chart(
        data,
        diffuse_cloud_label=diffuse_cloud_label,
        order=order,
    )
    figure.savefig(
        output_path,
        format=file_format,
        transparent=True,
        bbox_inches="tight",
        pad_inches=pad_inches,
    )
    _pyplot().close(figure)
    return output_path


def plot_individual_source_pie_chart(
    data: IndividualSourceData,
    *,
    ax=None,
    figure_size: tuple[float, float] = FIGURE_SIZE,
    ring_width: float = 0.10,
    ring_gap: float = 0.02,
    label_text_size: int = 14,
    percent_text_size: int = 12,
    remainder_color: str = "#EEEEEE",
    order: str = "count",
):
    """Plot the individual first-detection source concentric ring chart."""
    plt = _pyplot()

    if ax is None:
        figure = plt.figure(num="Individual Source Pie Chart", figsize=figure_size)
        ax = figure.add_subplot(111)
    else:
        figure = ax.figure

    categories = _ring_categories_by_order(data.categories, order)

    for index, category in enumerate(categories):
        radius = 1 - index * (ring_width + ring_gap)
        ax.pie(
            [category.fraction, 1.0 - category.fraction],
            colors=[category.color, remainder_color],
            radius=radius,
            startangle=_type_pie_start_angle(category.fraction),
            wedgeprops={
                "width": ring_width,
                "edgecolor": "white",
                "linewidth": 1,
            },
        )

    label_y_positions = (0.110, 0.160, 0.205, 0.255, 0.305)
    percent_y_positions = (0.870, 0.825, 0.775, 0.725, 0.680)
    for category, label_y, percent_y in zip(
        categories,
        label_y_positions,
        percent_y_positions,
    ):
        ax.annotate(
            category.label,
            xy=(0.5, label_y),
            xycoords="axes fraction",
            color=category.color,
            ha="center",
            size=label_text_size,
            weight="bold",
        )
        ax.annotate(
            f"{category.percent:.1f}%",
            xy=(0.51, percent_y),
            xycoords="axes fraction",
            color="white",
            ha="center",
            size=percent_text_size,
            weight="bold",
        )

    ax.set_aspect("equal")
    ax.axis("off")
    figure.tight_layout()
    return figure, ax


def write_individual_source_pie_chart(
    data: IndividualSourceData,
    output_path: str | Path,
    *,
    file_format: str | None = None,
    pad_inches: float = -0.65,
    order: str = "count",
) -> Path:
    """Write the individual-source ring chart and return the path."""
    output_path = Path(output_path)
    if file_format is None:
        file_format = output_path.suffix.lstrip(".") or "pdf"

    figure, _ = plot_individual_source_pie_chart(
        data,
        order=order,
    )
    figure.savefig(
        output_path,
        format=file_format,
        transparent=True,
        bbox_inches="tight",
        pad_inches=pad_inches,
    )
    _pyplot().close(figure)
    return output_path


def plot_molecule_type_by_source_type(
    data: MoleculeTypeBySourceData,
    *,
    figure_size: tuple[float, float] = (15, 12),
    label_text_size: int = 26,
    title_text_size: int = 30,
    legend_text_size: int = 26,
    wedge_alpha: float = 0.5,
    source_label_overrides: dict[str, str] | None = None,
):
    """Plot molecule-type pie charts within generalized source categories."""
    plt = _pyplot()
    from matplotlib.patches import Patch

    source_label_overrides = source_label_overrides or {}
    figure, axes = plt.subplots(
        2,
        2,
        num="Molecule Type by Source Type",
        figsize=figure_size,
    )

    for ax, source_category in zip(axes.flat, data.categories):
        nonzero_type_specs = [
            type_spec
            for type_spec in MOLECULE_TYPE_BY_SOURCE_TYPE_SPECS
            if source_category.counts_by_type[type_spec["key"]] > 0
        ]
        values = [
            source_category.counts_by_type[type_spec["key"]]
            for type_spec in nonzero_type_specs
        ]
        colors = [
            type_spec["color"]
            for type_spec in nonzero_type_specs
        ]

        if values:
            ax.pie(
                values,
                labels=[str(value) for value in values],
                colors=colors,
                labeldistance=0.8,
                textprops={"fontsize": label_text_size},
                wedgeprops={
                    "linewidth": 1.0,
                    "edgecolor": "black",
                    "alpha": wedge_alpha,
                },
            )
        title = source_label_overrides.get(
            source_category.key,
            source_category.label,
        )
        ax.annotate(
            title,
            xy=(0.5, 1.0),
            xycoords="axes fraction",
            color="black",
            ha="center",
            va="top",
            size=title_text_size,
            weight="bold",
        )
        ax.set_aspect("equal")

    legend_handles = [
        Patch(
            facecolor=type_spec["color"],
            edgecolor="black",
            alpha=wedge_alpha,
            label=type_spec["label"],
        )
        for type_spec in MOLECULE_TYPE_BY_SOURCE_TYPE_SPECS
    ]
    axes[0, 1].legend(
        handles=legend_handles,
        title="Molecule Types",
        loc="center left",
        bbox_to_anchor=(1, 0, 0.5, 1),
        fontsize=legend_text_size,
        title_fontsize=legend_text_size,
    )

    figure.tight_layout()
    figure.subplots_adjust(wspace=-0.3, hspace=0)
    return figure, axes


def write_molecule_type_by_source_type(
    data: MoleculeTypeBySourceData,
    output_path: str | Path,
    *,
    file_format: str | None = None,
    source_label_overrides: dict[str, str] | None = None,
) -> Path:
    """Write the molecule-type-by-source pie grid and return the path."""
    output_path = Path(output_path)
    if file_format is None:
        file_format = output_path.suffix.lstrip(".") or "pdf"

    figure, _ = plot_molecule_type_by_source_type(
        data,
        source_label_overrides=source_label_overrides,
    )
    figure.savefig(
        output_path,
        format=file_format,
        transparent=True,
        bbox_inches="tight",
    )
    _pyplot().close(figure)
    return output_path


def _molecule_type_source_display_label(
    category: MoleculeTypesForSource,
    overrides: dict[str, str],
) -> str:
    """Return the compact source label used in matrix-style figures."""
    if category.key in overrides:
        return overrides[category.key]
    labels = {
        "carbon_star": "Carbon Star",
        "dark_cloud": "Dark Cloud",
        "diffuse_cloud": "Diffuse Cloud",
        "sfr": "SFR",
    }
    return labels.get(category.key, category.label)


def _molecule_type_source_enrichment_factor(
    data: MoleculeTypeBySourceData,
    category: MoleculeTypesForSource,
    type_key: str,
) -> float:
    """Return fractional enrichment relative to the overall ISM/CSM mix."""
    if category.source_count == 0:
        return 0.0
    overall_count = data.overall_type_counts[type_key]
    if overall_count == 0 or data.molecule_count == 0:
        return 0.0
    source_fraction = category.counts_by_type[type_key] / category.source_count
    overall_fraction = overall_count / data.molecule_count
    return source_fraction / overall_fraction


def _format_enrichment_factor(factor: float) -> str:
    """Return a compact display label for a fractional enrichment factor."""
    if factor == 0:
        return "0.0x"
    if factor < 0.095:
        return "<0.1x"
    return f"{factor:.1f}x"


def plot_molecule_type_by_source_enrichment_matrix(
    data: MoleculeTypeBySourceData,
    *,
    figure_size: tuple[float, float] = (12, 8),
    axes_bounds: tuple[float, float, float, float] = (0.205, 0.15, 0.62, 0.775),
    colorbar_bounds: tuple[float, float, float, float] = (0.86, 0.15, 0.035, 0.775),
    source_label_overrides: dict[str, str] | None = None,
    text_size: int = FIGURE_TEXT_SIZE,
    factor_text_size: int = 21,
    count_text_size: int = 17,
    depletion_limit: float = 0.25,
    enrichment_limit: float = 4.0,
):
    """Plot the production molecule-type/source enrichment matrix.

    Cell values are fractional enrichments relative to the overall secure
    ISM/CSM molecule-type mix. Color is log-scaled and centered on ``1x``.
    """
    plt = _pyplot()
    from matplotlib import colors
    from matplotlib.colors import LinearSegmentedColormap

    if depletion_limit <= 0 or enrichment_limit <= 1:
        raise ValueError("Use depletion_limit > 0 and enrichment_limit > 1.")

    source_label_overrides = source_label_overrides or {}
    type_specs_by_key = {
        type_spec["key"]: type_spec
        for type_spec in MOLECULE_TYPE_BY_SOURCE_TYPE_SPECS
    }
    type_keys = list(MOLECULE_TYPE_BY_SOURCE_MATRIX_TYPE_KEYS)
    type_labels = [
        type_specs_by_key[type_key]["label"]
        for type_key in type_keys
    ]
    source_labels = [
        _molecule_type_source_display_label(category, source_label_overrides)
        for category in data.categories
    ]

    enrichment_factors = np.array(
        [
            [
                _molecule_type_source_enrichment_factor(data, category, type_key)
                for type_key in type_keys
            ]
            for category in data.categories
        ],
        dtype=float,
    )

    min_log = np.log2(depletion_limit)
    max_log = np.log2(enrichment_limit)
    log_enrichment = np.full_like(enrichment_factors, min_log, dtype=float)
    positive_mask = enrichment_factors > 0
    log_enrichment[positive_mask] = np.log2(enrichment_factors[positive_mask])
    log_enrichment = np.clip(log_enrichment, min_log, max_log)

    colormap = LinearSegmentedColormap.from_list(
        "astromol_gray_white_blue",
        MOLECULE_TYPE_BY_SOURCE_ENRICHMENT_COLORS,
    )
    norm = colors.TwoSlopeNorm(vmin=min_log, vcenter=0.0, vmax=max_log)

    figure = plt.figure(
        num="Molecule Type by Source Enrichment Matrix",
        figsize=figure_size,
    )
    ax = figure.add_axes(axes_bounds)
    colorbar_ax = figure.add_axes(colorbar_bounds)

    image = ax.imshow(log_enrichment, cmap=colormap, norm=norm, aspect="equal")
    ax.set_xticks(np.arange(len(type_labels)))
    ax.set_xticklabels(type_labels, fontsize=text_size)
    ax.set_yticks(np.arange(len(source_labels)))
    ax.set_yticklabels(source_labels, fontsize=text_size)
    ax.set_xlabel("Molecule Type", fontsize=text_size)
    ax.set_ylabel("First-Detection Source Type", fontsize=text_size)
    ax.tick_params(
        axis="both",
        which="both",
        direction="in",
        length=8,
        width=FIGURE_TICK_WIDTH,
        colors="black",
    )
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")

    for row, category in enumerate(data.categories):
        for col, type_key in enumerate(type_keys):
            count = category.counts_by_type[type_key]
            factor = enrichment_factors[row, col]
            ax.text(
                col,
                row - 0.055,
                _format_enrichment_factor(factor),
                ha="center",
                va="center",
                fontsize=factor_text_size,
                color="black",
                fontstyle="italic",
                weight="semibold",
            )
            ax.text(
                col,
                row + 0.155,
                f"n={count}",
                ha="center",
                va="center",
                fontsize=count_text_size,
                color="black",
                alpha=0.82,
            )

    ax.set_xticks(np.arange(-0.5, len(type_labels), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(source_labels), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.1)
    ax.tick_params(which="minor", bottom=False, left=False)

    colorbar = figure.colorbar(image, cax=colorbar_ax)
    colorbar.set_label(
        "Fractional Enrichment vs Overall ISM/CSM",
        fontsize=18,
    )
    colorbar.ax.tick_params(labelsize=16, colors="black")
    colorbar.set_ticks([min_log, -1, 0, 1, max_log])
    colorbar.set_ticklabels(
        [
            f"{depletion_limit:g}x",
            "0.5x",
            "1x",
            "2x",
            f"{enrichment_limit:g}x",
        ]
    )
    return figure, ax


def write_molecule_type_by_source_enrichment_matrix(
    data: MoleculeTypeBySourceData,
    output_path: str | Path,
    *,
    file_format: str | None = None,
    source_label_overrides: dict[str, str] | None = None,
) -> Path:
    """Write the production molecule-type/source enrichment matrix."""
    output_path = Path(output_path)
    if file_format is None:
        file_format = output_path.suffix.lstrip(".") or "pdf"

    figure, _ = plot_molecule_type_by_source_enrichment_matrix(
        data,
        source_label_overrides=source_label_overrides,
    )
    figure.savefig(
        output_path,
        format=file_format,
        transparent=True,
        bbox_inches="tight",
    )
    _pyplot().close(figure)
    return output_path


def _kde_values(
    masses: tuple[float, ...],
    x_values: np.ndarray,
    *,
    bandwidth: float,
) -> np.ndarray:
    """Return Gaussian KDE values using the legacy bandwidth convention."""
    from scipy.stats import gaussian_kde

    density = gaussian_kde(masses)
    density.covariance_factor = lambda: bandwidth
    density._compute_covariance()
    return density(x_values)


def _underline_axes_text(
    ax,
    text_artist,
    *,
    y_offset: float = 0.008,
    linewidth: float = 1.0,
) -> None:
    """Draw an underline below axes-fraction text using its rendered width."""
    figure = ax.figure
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    text_bbox = text_artist.get_window_extent(renderer=renderer)
    axes_inverse = ax.transAxes.inverted()
    x0, y0 = axes_inverse.transform((text_bbox.x0, text_bbox.y0))
    x1, _ = axes_inverse.transform((text_bbox.x1, text_bbox.y0))
    ax.plot(
        [x0, x1],
        [y0 - y_offset, y0 - y_offset],
        transform=ax.transAxes,
        color="black",
        linewidth=linewidth,
        clip_on=False,
    )


def _category_by_key(categories) -> dict[str, object]:
    """Return categories keyed by their short identifier."""
    return {
        category.key: category
        for category in categories
    }


def _auto_du_source_count_position(
    ax,
    curve: _DUBySourceTypeCurve,
    curves: tuple[_DUBySourceTypeCurve, ...],
    placed_boxes: list[tuple[float, float, float, float]],
    *,
    text_size: int,
    legend_box_axes: tuple[float, float, float, float],
) -> tuple[float, float]:
    """Choose a readable count-label position near a DU/source KDE curve."""
    label = str(curve.category.count)
    label_width, label_height = _mass_label_size_data(ax, label, text_size)
    x_min, x_max = ax.get_xlim()
    y_min, y_max = ax.get_ylim()
    x_span = x_max - x_min
    y_span = y_max - y_min
    x_pad = x_span * 0.018
    y_pad = y_span * 0.025
    y_margin = y_span * 0.004

    peak_x = curve.peak_x
    peak_y = curve.peak_y
    top_label_y = y_max - label_height - y_pad
    shoulder_x = min(peak_x + x_span * 0.10, x_max - label_width - x_pad)
    shoulder_y = float(np.interp(shoulder_x, curve.x_values, curve.density))
    tail_x = min(peak_x + x_span * 0.38, x_max - label_width - x_pad)
    tail_y = float(np.interp(tail_x, curve.x_values, curve.density))
    far_tail_x = min(peak_x + x_span * 0.52, x_max - label_width - x_pad)
    far_tail_y = float(np.interp(far_tail_x, curve.x_values, curve.density))

    candidate_positions = [
        (peak_x + x_pad, peak_y - label_height * 0.35),
        (peak_x + x_pad, peak_y + y_pad),
        (shoulder_x, shoulder_y + y_pad),
        (shoulder_x, shoulder_y - label_height * 0.35),
        (peak_x - label_width - x_pad, peak_y - label_height * 0.35),
        (peak_x - label_width - x_pad, peak_y + y_pad),
        (tail_x, tail_y + y_pad),
        (tail_x, tail_y + y_pad * 2.5),
        (far_tail_x, far_tail_y + y_pad),
        (far_tail_x, far_tail_y + y_pad * 2.5),
    ]
    if peak_x < x_min + x_span * 0.08:
        candidate_positions = [
            (peak_x + x_pad, top_label_y),
            (peak_x + x_pad * 2.5, top_label_y),
            *candidate_positions,
        ]

    best_position = None
    best_score = float("inf")
    for x0, y0 in candidate_positions:
        box = (x0, y0, x0 + label_width, y0 + label_height)
        score = (
            abs((x0 + label_width / 2) - peak_x) / max(x_span, 1)
            + abs((y0 + label_height / 2) - peak_y) / max(y_span, 1e-12)
        )
        if box[0] < x_min or box[2] > x_max or box[1] < y_min or box[3] > y_max:
            score += 100

        axes_box = _data_box_to_axes_box(ax, box)
        if _boxes_overlap(axes_box, legend_box_axes):
            score += 50
        for placed_box in placed_boxes:
            if _boxes_overlap(box, placed_box):
                score += 30
        for other_curve in curves:
            if _curve_line_overlaps_data_box(other_curve, box, y_margin=y_margin):
                score += 15
            if (
                other_curve.category.key != curve.category.key
                and _curve_fill_overlaps_data_box(other_curve, box, y_margin=y_margin)
            ):
                score += 25
        if score < best_score:
            best_score = score
            best_position = (x0, y0)

    if best_position is None:
        best_position = (
            min(max(peak_x + x_pad, x_min), x_max - label_width),
            min(max(peak_y + y_pad, y_min), y_max - label_height),
        )
    return best_position


def _du_source_count_positions(
    ax,
    curves: tuple[_DUBySourceTypeCurve, ...],
    *,
    text_size: int,
    count_label_overrides: dict[str, tuple[float, float]] | None,
    legend_box_axes: tuple[float, float, float, float],
) -> dict[str, tuple[float, float]]:
    """Return count-label positions keyed by source category."""
    positions: dict[str, tuple[float, float]] = {}
    placed_boxes: list[tuple[float, float, float, float]] = []
    for curve in sorted(curves, key=lambda item: item.peak_y, reverse=True):
        if count_label_overrides and curve.category.key in count_label_overrides:
            position = count_label_overrides[curve.category.key]
        else:
            position = _auto_du_source_count_position(
                ax,
                curve,
                curves,
                placed_boxes,
                text_size=text_size,
                legend_box_axes=legend_box_axes,
            )

        label_width, label_height = _mass_label_size_data(
            ax,
            str(curve.category.count),
            text_size,
        )
        placed_boxes.append(
            (
                position[0],
                position[1],
                position[0] + label_width,
                position[1] + label_height,
            )
        )
        positions[curve.category.key] = position
    return positions


def plot_du_by_source_type(
    data: DUBySourceTypeData,
    *,
    ax=None,
    bandwidth: float = 0.5,
    figure_size: tuple[float, float] = FIGURE_SIZE,
    text_size: int = FIGURE_TEXT_SIZE,
    source_label_overrides: dict[str, str] | None = None,
    count_label_overrides: dict[str, tuple[float, float]] | None = None,
    ylabel: str = "Probability Density Estimate",
):
    """Plot DU KDE distributions by generalized first-detection source type."""
    plt = _pyplot()
    source_label_overrides = source_label_overrides or {}

    if ax is None:
        figure = plt.figure(num="DU by Source Type", figsize=figure_size)
        ax = figure.add_subplot(111)
    else:
        figure = ax.figure

    x_values = np.arange(0, 15, 0.1)
    categories_by_key = _category_by_key(data.categories)
    max_density = 0.0
    curves = []
    for source_spec in DU_BY_SOURCE_TYPE_SPECS:
        category = categories_by_key[source_spec["key"]]
        if category.count < 2 or len(set(category.values)) < 2:
            continue

        density = _kde_values(
            category.values,
            x_values,
            bandwidth=bandwidth,
        )
        max_density = max(max_density, float(np.max(density)))
        curves.append(
            _DUBySourceTypeCurve(
                category=category,
                x_values=x_values,
                density=density,
            )
        )
        ax.plot(x_values, density, color=category.color)
        ax.fill_between(
            x_values,
            density,
            0,
            facecolor=category.color,
            alpha=0.25,
            zorder=4,
        )

    header = ax.annotate(
        "Source Types",
        xy=(0.97, 0.96),
        xycoords="axes fraction",
        color="black",
        ha="right",
        va="top",
        size=text_size,
    )
    for row, source_key in enumerate(DU_BY_SOURCE_TYPE_LEGEND_ORDER, start=1):
        category = categories_by_key[source_key]
        label = source_label_overrides.get(category.key, category.label)
        ax.annotate(
            label,
            xy=(0.97, 0.96 - 0.06 * row),
            xycoords="axes fraction",
            color=category.color,
            ha="right",
            va="top",
            size=text_size,
        )

    ax.set_xlabel("Degree of Unsaturation", fontsize=text_size)
    ax.set_ylabel(ylabel, fontsize=text_size)
    ax.set_xlim(0, 15)
    ax.set_ylim(0, max(0.4, float(np.ceil(max_density * 10) / 10)))
    ax.tick_params(
        axis="both",
        which="both",
        direction="in",
        length=5,
        width=FIGURE_TICK_WIDTH,
        labelsize=text_size,
        colors="black",
    )
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
    ax.grid(False)
    ax.set_title("")

    figure.tight_layout()
    count_positions = _du_source_count_positions(
        ax,
        tuple(curves),
        text_size=text_size,
        count_label_overrides=count_label_overrides,
        legend_box_axes=(0.56, 0.70, 1.0, 1.0),
    )
    for source_key in DU_BY_SOURCE_TYPE_LEGEND_ORDER:
        category = categories_by_key[source_key]
        if category.key not in count_positions:
            continue
        ax.annotate(
            str(category.count),
            xy=count_positions[category.key],
            xycoords="data",
            ha="left",
            va="bottom",
            color=category.color,
            size=text_size,
        )
    _underline_axes_text(ax, header)
    return figure, ax


def write_du_by_source_type(
    data: DUBySourceTypeData,
    output_path: str | Path,
    *,
    bandwidth: float = 0.5,
    source_label_overrides: dict[str, str] | None = None,
    count_label_overrides: dict[str, tuple[float, float]] | None = None,
    ylabel: str = "Probability Density Estimate",
    file_format: str | None = None,
) -> Path:
    """Write the DU-by-source KDE plot and return the output path."""
    output_path = Path(output_path)
    if file_format is None:
        file_format = output_path.suffix.lstrip(".") or "pdf"

    figure, _ = plot_du_by_source_type(
        data,
        bandwidth=bandwidth,
        source_label_overrides=source_label_overrides,
        count_label_overrides=count_label_overrides,
        ylabel=ylabel,
    )
    figure.savefig(
        output_path,
        format=file_format,
        transparent=True,
        bbox_inches="tight",
        pad_inches=0,
    )
    _pyplot().close(figure)
    return output_path


def plot_relative_du_by_source_type(
    data: RelativeDUBySourceTypeData,
    *,
    axes=None,
    bandwidth: float = 0.5,
    figure_size: tuple[float, float] = FIGURE_SIZE,
    text_size: int = 18,
    label_size: int = FIGURE_TEXT_SIZE,
    source_label_overrides: dict[str, str] | None = None,
):
    """Plot the legacy relative-DU KDE panel by source type."""
    plt = _pyplot()
    source_label_overrides = source_label_overrides or {}

    if axes is None:
        figure, axes_array = plt.subplots(
            2,
            2,
            num="Relative DU by Source Type",
            figsize=figure_size,
        )
        axes_list = tuple(axes_array.flat)
    else:
        axes_list = tuple(axes)
        figure = axes_list[0].figure

    categories_by_key = _category_by_key(data.categories)
    x_values = np.arange(0, 1, 0.01)
    for index, (ax, source_key) in enumerate(
        zip(axes_list, RELATIVE_DU_BY_SOURCE_TYPE_PANEL_ORDER)
    ):
        category = categories_by_key[source_key]
        if category.count >= 2 and len(set(category.values)) >= 2:
            density = _kde_values(
                category.values,
                x_values,
                bandwidth=bandwidth,
            )
            ax.plot(x_values, density, color=category.color)
            ax.fill_between(
                x_values,
                density,
                0,
                facecolor=category.color,
                alpha=0.25,
                zorder=4,
            )

        ax.annotate(
            source_label_overrides.get(category.key, category.label),
            xy=(0.04, 0.96),
            xycoords="axes fraction",
            ha="left",
            va="top",
            size=label_size,
            color=category.color,
        )
        ax.set_xlim(0, 1)
        ax.tick_params(
            axis="both",
            which="both",
            direction="in",
            length=5,
            width=FIGURE_TICK_WIDTH,
            labelsize=text_size,
            colors="black",
        )
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_color("black")
        if index != 2:
            ax.xaxis.set_visible(False)
            ax.yaxis.set_visible(False)

    axes_list[2].set_xlabel("Relative Degree of Unsaturation", fontsize=text_size)
    axes_list[2].set_ylabel("Probability Density Estimate", fontsize=text_size)
    figure.subplots_adjust(wspace=0, hspace=0)
    figure.tight_layout()
    return figure, axes_list


def write_relative_du_by_source_type(
    data: RelativeDUBySourceTypeData,
    output_path: str | Path,
    *,
    bandwidth: float = 0.5,
    source_label_overrides: dict[str, str] | None = None,
    file_format: str | None = None,
) -> Path:
    """Write the relative-DU-by-source KDE plot and return the output path."""
    output_path = Path(output_path)
    if file_format is None:
        file_format = output_path.suffix.lstrip(".") or "pdf"

    figure, _ = plot_relative_du_by_source_type(
        data,
        bandwidth=bandwidth,
        source_label_overrides=source_label_overrides,
    )
    figure.savefig(
        output_path,
        format=file_format,
        transparent=True,
        bbox_inches="tight",
    )
    _pyplot().close(figure)
    return output_path


def plot_mass_by_source_type(
    data: MassBySourceTypeData,
    *,
    ax=None,
    bandwidth: float = 0.5,
    figure_size: tuple[float, float] = FIGURE_SIZE,
    text_size: int = FIGURE_TEXT_SIZE,
    source_label_overrides: dict[str, str] | None = None,
):
    """Plot legacy molecular-mass KDE distributions by source type."""
    plt = _pyplot()
    source_label_overrides = source_label_overrides or {}

    if ax is None:
        figure = plt.figure(num="Mass by Source Type", figsize=figure_size)
        ax = figure.add_subplot(111)
    else:
        figure = ax.figure

    categories_by_key = _category_by_key(data.categories)
    x_min, x_max = data.mass_range
    x_values = np.arange(0, int(np.ceil(x_max)) + 1, 1)
    max_density = 0.0
    for source_spec in DU_BY_SOURCE_TYPE_SPECS:
        category = categories_by_key[source_spec["key"]]
        if category.count < 2 or len(set(category.masses)) < 2:
            continue

        density = _kde_values(
            category.masses,
            x_values,
            bandwidth=bandwidth,
        )
        max_density = max(max_density, float(np.max(density)))
        ax.plot(x_values, density, color=category.color)
        ax.fill_between(
            x_values,
            density,
            0,
            facecolor=category.color,
            alpha=0.25,
            zorder=4,
        )

    header = ax.annotate(
        "Source Types",
        xy=(0.97, 0.96),
        xycoords="axes fraction",
        ha="right",
        va="top",
        size=text_size,
        color="black",
    )
    for row, source_key in enumerate(DU_BY_SOURCE_TYPE_LEGEND_ORDER, start=1):
        category = categories_by_key[source_key]
        ax.annotate(
            source_label_overrides.get(category.key, category.label),
            xy=(0.97, 0.96 - 0.06 * row),
            xycoords="axes fraction",
            color=category.color,
            ha="right",
            va="top",
            size=text_size,
        )

    ax.set_xlabel("Molecular Mass (amu)", fontsize=text_size)
    ax.set_ylabel("Probability Density Estimate", fontsize=text_size)
    y_top = max(0.03, float(np.ceil(max_density * 1000) / 1000))
    ax.set_xlim(x_min, x_max)
    ax.set_xticks(np.arange(25, int(np.ceil(x_max)) + 1, 25))
    ax.set_ylim(0, y_top)
    ax.set_yticks(np.arange(0, y_top, 0.005))
    ax.tick_params(
        axis="both",
        which="both",
        direction="in",
        length=5,
        width=FIGURE_TICK_WIDTH,
        labelsize=text_size,
        colors="black",
    )
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
    ax.grid(False)
    ax.set_title("")

    _underline_axes_text(ax, header)
    return figure, ax


def write_mass_by_source_type(
    data: MassBySourceTypeData,
    output_path: str | Path,
    *,
    bandwidth: float = 0.5,
    source_label_overrides: dict[str, str] | None = None,
    file_format: str | None = None,
) -> Path:
    """Write the mass-by-source KDE plot and return the output path."""
    output_path = Path(output_path)
    if file_format is None:
        file_format = output_path.suffix.lstrip(".") or "pdf"

    figure, _ = plot_mass_by_source_type(
        data,
        bandwidth=bandwidth,
        source_label_overrides=source_label_overrides,
    )
    figure.savefig(
        output_path,
        format=file_format,
        transparent=True,
        bbox_inches="tight",
        pad_inches=0,
    )
    _pyplot().close(figure)
    return output_path


def _wavelength_display_label(wavelength: str) -> str:
    """Return manuscript label for a wavelength category."""
    if wavelength == "UV-Vis":
        return "UV/Vis"
    return wavelength


def _pie_percentage_labels(values: list[int]) -> list[str]:
    """Return one-decimal percentage labels for non-empty pie values."""
    total = sum(values)
    if total == 0:
        return []
    return [
        f"{100 * value / total:.1f}%"
        for value in values
    ]


def plot_wavelength_by_source_type(
    data: WavelengthBySourceTypeData,
    *,
    axes=None,
    figure_size: tuple[float, float] = (15, 12),
    text_size: int = 26,
    panel_label_size: int = 30,
    source_label_overrides: dict[str, str] | None = None,
):
    """Plot the legacy wavelength-by-source-type pie grid."""
    plt = _pyplot()
    source_label_overrides = source_label_overrides or {}

    if axes is None:
        figure, axes_array = plt.subplots(
            2,
            2,
            num="Wavelength by Source Type",
            figsize=figure_size,
        )
        axes_list = tuple(axes_array.flat)
    else:
        axes_list = tuple(axes)
        figure = axes_list[0].figure

    categories_by_key = _category_by_key(data.categories)
    for ax, source_key in zip(axes_list, WAVES_BY_SOURCE_TYPE_PANEL_ORDER):
        category = categories_by_key[source_key]
        display_wavelengths = WAVES_BY_SOURCE_TYPE_DISPLAY_WAVELENGTHS[source_key]
        pie_values = []
        pie_colors = []
        for wavelength in display_wavelengths:
            count = category.count_for_display_wavelength(wavelength)
            if count <= 0:
                continue
            pie_values.append(count)
            pie_colors.append(WAVES_BY_SOURCE_TYPE_COLORS[wavelength])

        if pie_values:
            ax.pie(
                pie_values,
                labels=_pie_percentage_labels(pie_values),
                colors=pie_colors,
                labeldistance=1.1,
                textprops={"fontsize": text_size, "color": "black"},
                wedgeprops={
                    "linewidth": 1.0,
                    "edgecolor": "black",
                    "alpha": 0.5,
                },
            )
        panel_label = source_label_overrides.get(
            source_key,
            WAVES_BY_SOURCE_TYPE_PANEL_LABELS[source_key],
        )
        ax.annotate(
            panel_label,
            xy=(0.5, 1.05),
            xycoords="axes fraction",
            color="black",
            ha="center",
            va="top",
            size=panel_label_size,
            fontweight="bold",
        )
        ax.set_aspect("equal")
        ax.set_title("")

    from matplotlib.patches import Patch

    legend_wavelengths = ("cm", "mm", "sub-mm", "IR", "UV-Vis")
    legend_handles = [
        Patch(
            facecolor=WAVES_BY_SOURCE_TYPE_COLORS[wavelength],
            edgecolor="black",
            linewidth=1.0,
            alpha=0.5,
        )
        for wavelength in legend_wavelengths
    ]
    legend = axes_list[1].legend(
        legend_handles,
        [_wavelength_display_label(wavelength) for wavelength in legend_wavelengths],
        title="Wavelengths",
        loc="center left",
        bbox_to_anchor=(1, 0, 0.5, 1),
        frameon=True,
        fontsize=text_size,
        title_fontsize=text_size,
    )
    legend.get_frame().set_edgecolor("#cccccc")

    figure.tight_layout()
    figure.subplots_adjust(wspace=0.1, hspace=0.1)
    return figure, axes_list


def write_wavelength_by_source_type(
    data: WavelengthBySourceTypeData,
    output_path: str | Path,
    *,
    source_label_overrides: dict[str, str] | None = None,
    file_format: str | None = None,
) -> Path:
    """Write the wavelength-by-source-type pie grid and return the path."""
    output_path = Path(output_path)
    if file_format is None:
        file_format = output_path.suffix.lstrip(".") or "pdf"

    figure, _ = plot_wavelength_by_source_type(
        data,
        source_label_overrides=source_label_overrides,
    )
    figure.savefig(
        output_path,
        format=file_format,
        transparent=True,
        bbox_inches="tight",
    )
    _pyplot().close(figure)
    return output_path


def plot_wavelength_by_source_type_stacked_bar(
    data: WavelengthBySourceTypeData,
    *,
    ax=None,
    figure_size: tuple[float, float] = (8.0, 4.8),
    text_size: float = 18.0,
    tick_size: float = 16.0,
    segment_text_size: float = 13.0,
    show_segment_labels: bool = False,
    source_label_overrides: dict[str, str] | None = None,
):
    """Plot wavelength shares by source type as normalized stacked bars."""
    plt = _pyplot()
    source_label_overrides = source_label_overrides or {}

    if ax is None:
        figure = plt.figure(
            num="Wavelength Shares by Source Type",
            figsize=figure_size,
        )
        ax = figure.add_subplot(111)
    else:
        figure = ax.figure

    categories_by_key = _category_by_key(data.categories)
    categories = tuple(
        categories_by_key[key]
        for key in DU_BY_SOURCE_TYPE_BOXPLOT_ORDER
        if key in categories_by_key
    )
    display_wavelengths = ("cm", "mm", "sub-mm", "IR", "UV-Vis")
    counts_matrix = np.array(
        [
            [
                category.count_for_display_wavelength(wavelength)
                for wavelength in display_wavelengths
            ]
            for category in categories
        ],
        dtype=float,
    )
    totals = counts_matrix.sum(axis=1)
    percent_matrix = np.divide(
        counts_matrix,
        totals[:, np.newaxis],
        out=np.zeros_like(counts_matrix),
        where=totals[:, np.newaxis] != 0,
    ) * 100

    y_positions = np.arange(len(categories), 0, -1)
    left = np.zeros(len(categories), dtype=float)
    for col, wavelength in enumerate(display_wavelengths):
        percentages = percent_matrix[:, col]
        bars = ax.barh(
            y_positions,
            percentages,
            left=left,
            height=0.72,
            color=WAVES_BY_SOURCE_TYPE_COLORS[wavelength],
            edgecolor="black",
            linewidth=0.8,
            alpha=0.82,
            label=_wavelength_display_label(wavelength),
            zorder=3,
        )
        if show_segment_labels:
            for row, (bar, percentage) in enumerate(zip(bars, percentages)):
                if percentage < 7.0:
                    continue
                text_color = "white" if wavelength == "IR" else "black"
                ax.text(
                    left[row] + percentage / 2,
                    bar.get_y() + bar.get_height() / 2,
                    f"{percentage:.0f}%",
                    ha="center",
                    va="center",
                    fontsize=segment_text_size,
                    color=text_color,
                    zorder=4,
                )
        left += percentages

    from matplotlib.transforms import blended_transform_factory

    count_label_transform = blended_transform_factory(ax.transAxes, ax.transData)
    for y_position, total in zip(y_positions, totals.astype(int)):
        ax.text(
            1.015,
            y_position,
            f"n={total}",
            ha="left",
            va="center",
            fontsize=tick_size,
            color="#444444",
            transform=count_label_transform,
            clip_on=False,
        )

    ax.set_yticks(y_positions)
    ax.set_yticklabels(
        [
            source_label_overrides.get(category.key, category.label)
            for category in categories
        ],
        fontsize=tick_size,
    )
    ax.set_xlabel("First-detection wavelength credits (%)", fontsize=text_size)
    ax.set_xlim(0, 100)
    ax.set_ylim(0.45, len(categories) + 0.55)
    ax.set_xticks(np.arange(0, 101, 20))
    ax.tick_params(
        axis="both",
        which="both",
        direction="in",
        length=5,
        width=FIGURE_TICK_WIDTH,
        labelsize=tick_size,
        colors="black",
    )
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
    ax.grid(False)
    ax.set_title("")
    ax.legend(
        loc="lower center",
        bbox_to_anchor=(0.5, 1.02),
        ncol=len(display_wavelengths),
        frameon=False,
        fontsize=tick_size,
        handlelength=1.6,
        columnspacing=1.2,
    )
    figure.tight_layout()
    return figure, ax


def write_wavelength_by_source_type_stacked_bar(
    data: WavelengthBySourceTypeData,
    output_path: str | Path,
    *,
    source_label_overrides: dict[str, str] | None = None,
    file_format: str | None = None,
) -> Path:
    """Write the production wavelength/source stacked bar chart."""
    output_path = Path(output_path)
    if file_format is None:
        file_format = output_path.suffix.lstrip(".") or "pdf"

    figure, _ = plot_wavelength_by_source_type_stacked_bar(
        data,
        source_label_overrides=source_label_overrides,
    )
    figure.savefig(
        output_path,
        format=file_format,
        transparent=True,
        bbox_inches="tight",
    )
    _pyplot().close(figure)
    return output_path


def _data_box_to_axes_box(ax, box: tuple[float, float, float, float]) -> tuple[float, float, float, float]:
    """Convert a data-coordinate box to axes-fraction coordinates."""
    x0, y0, x1, y1 = box
    axes_inverse = ax.transAxes.inverted()
    ax_x0, ax_y0 = axes_inverse.transform(ax.transData.transform((x0, y0)))
    ax_x1, ax_y1 = axes_inverse.transform(ax.transData.transform((x1, y1)))
    return ax_x0, ax_y0, ax_x1, ax_y1


def _boxes_overlap(
    first: tuple[float, float, float, float],
    second: tuple[float, float, float, float],
) -> bool:
    """Return True when two ``(x0, y0, x1, y1)`` boxes overlap."""
    return not (
        first[2] <= second[0]
        or second[2] <= first[0]
        or first[3] <= second[1]
        or second[3] <= first[1]
    )


def _mass_label_size_data(ax, label: str, text_size: int) -> tuple[float, float]:
    """Estimate a count-label bounding box in data coordinates."""
    figure = ax.figure
    figure.canvas.draw()
    axes_bbox = ax.get_window_extent()
    x_min, x_max = ax.get_xlim()
    y_min, y_max = ax.get_ylim()
    font_px = text_size * figure.dpi / 72.0
    label_width_px = max(len(label), 1) * font_px * 0.62
    label_height_px = font_px * 1.12
    label_width_data = (x_max - x_min) * label_width_px / axes_bbox.width
    label_height_data = (y_max - y_min) * label_height_px / axes_bbox.height
    return label_width_data, label_height_data


def _curve_line_overlaps_data_box(
    curve: _MassByWavelengthCurve,
    box: tuple[float, float, float, float],
    *,
    y_margin: float,
) -> bool:
    """Return True if a curve line passes through a label candidate box."""
    x0, y0, x1, y1 = box
    if x1 < curve.x_values[0] or x0 > curve.x_values[-1]:
        return False
    sample_x = np.linspace(max(x0, curve.x_values[0]), min(x1, curve.x_values[-1]), 25)
    sample_y = np.interp(sample_x, curve.x_values, curve.density)
    return bool(np.any((sample_y >= y0 - y_margin) & (sample_y <= y1 + y_margin)))


def _curve_fill_overlaps_data_box(
    curve: _MassByWavelengthCurve,
    box: tuple[float, float, float, float],
    *,
    y_margin: float,
) -> bool:
    """Return True if a label candidate sits inside another filled KDE area."""
    x0, y0, x1, _ = box
    if x1 < curve.x_values[0] or x0 > curve.x_values[-1]:
        return False
    sample_x = np.linspace(max(x0, curve.x_values[0]), min(x1, curve.x_values[-1]), 25)
    sample_y = np.interp(sample_x, curve.x_values, curve.density)
    return bool(np.any(sample_y >= y0 - y_margin))


def _mass_label_override(
    series: MassByWavelengthSeries,
    label_overrides: dict[str, tuple[float, float]] | None,
) -> tuple[float, float] | None:
    """Return a manual label position override for a wavelength series."""
    if not label_overrides:
        return None
    return label_overrides.get(series.key) or label_overrides.get(series.label)


def _auto_mass_label_position(
    ax,
    curve: _MassByWavelengthCurve,
    curves: tuple[_MassByWavelengthCurve, ...],
    placed_boxes: list[tuple[float, float, float, float]],
    *,
    text_size: int,
    legend_box_axes: tuple[float, float, float, float],
) -> tuple[float, float]:
    """Choose a near-peak label position that avoids other labels and curves."""
    label = str(curve.series.count)
    label_width, label_height = _mass_label_size_data(ax, label, text_size)
    x_min, x_max = ax.get_xlim()
    y_min, y_max = ax.get_ylim()
    x_span = x_max - x_min
    y_span = y_max - y_min
    x_pad = x_span * 0.015
    y_pad = y_span * 0.025
    y_margin = y_span * 0.004

    peak_x = curve.peak_x
    peak_y = curve.peak_y
    right_shoulder_x = min(peak_x + x_span * 0.06, x_max - label_width - x_pad)
    far_right_x = min(curve.x_values[-1] - label_width - x_pad, x_max - label_width - x_pad)
    right_shoulder_y = float(np.interp(right_shoulder_x, curve.x_values, curve.density))
    far_right_y = float(np.interp(far_right_x, curve.x_values, curve.density))
    low_flat_trace = peak_y < y_span * 0.12
    candidate_positions = [
        (peak_x + x_pad, peak_y + y_pad),
        (peak_x + x_pad, peak_y + y_pad * 2.8),
        (right_shoulder_x, right_shoulder_y + y_pad),
        (right_shoulder_x, right_shoulder_y + y_pad * 2.2),
        (peak_x - label_width - x_pad, peak_y + y_pad),
        (peak_x + x_pad, peak_y - label_height * 0.35),
        (peak_x - label_width - x_pad, peak_y - label_height * 0.35),
        (peak_x + x_pad * 2.5, peak_y - label_height * 1.2),
        (peak_x - label_width - x_pad * 2.5, peak_y - label_height * 1.2),
        (peak_x + x_pad * 4.5, peak_y + y_pad * 2.0),
        (peak_x - label_width - x_pad * 4.5, peak_y + y_pad * 2.0),
    ]
    if low_flat_trace:
        candidate_positions.extend(
            [
                (far_right_x, far_right_y + y_pad),
                (far_right_x, peak_y + y_pad),
            ]
        )

    best_position = None
    best_score = float("inf")
    preferred_x = peak_x
    if low_flat_trace:
        preferred_x = far_right_x
    for x0, y0 in candidate_positions:
        box = (x0, y0, x0 + label_width, y0 + label_height)
        score = (
            abs((x0 + label_width / 2) - preferred_x) / max(x_span, 1)
            + abs((y0 + label_height / 2) - peak_y) / max(y_span, 1e-12)
        )
        if box[0] < x_min or box[2] > x_max or box[1] < y_min or box[3] > y_max:
            score += 100
        axes_box = _data_box_to_axes_box(ax, box)
        if _boxes_overlap(axes_box, legend_box_axes):
            score += 50
        for placed_box in placed_boxes:
            if _boxes_overlap(box, placed_box):
                score += 25
        for other_curve in curves:
            if _curve_line_overlaps_data_box(other_curve, box, y_margin=y_margin):
                score += 18
            if (
                other_curve.series.key != curve.series.key
                and _curve_fill_overlaps_data_box(other_curve, box, y_margin=y_margin)
            ):
                score += 30
        if score < best_score:
            best_score = score
            best_position = (x0, y0)

    if best_position is None:
        clipped_x = min(max(peak_x + x_pad, x_min), x_max - label_width)
        clipped_y = min(max(peak_y + y_pad, y_min), y_max - label_height)
        best_position = (clipped_x, clipped_y)

    return best_position


def _mass_label_positions(
    ax,
    curves: tuple[_MassByWavelengthCurve, ...],
    *,
    text_size: int,
    label_mode: str,
    label_overrides: dict[str, tuple[float, float]] | None,
    legend_box_axes: tuple[float, float, float, float],
) -> dict[str, tuple[float, float]]:
    """Return count-label positions keyed by wavelength series key."""
    positions: dict[str, tuple[float, float]] = {}
    if label_mode not in {"fixed", "auto"}:
        raise ValueError("label_mode must be 'fixed' or 'auto'.")

    placed_boxes: list[tuple[float, float, float, float]] = []
    ordered_curves = sorted(curves, key=lambda item: item.peak_y, reverse=True)
    for curve in ordered_curves:
        override = _mass_label_override(curve.series, label_overrides)
        if override is not None:
            position = override
        elif label_mode == "fixed":
            position = curve.series.annotation_xy
        else:
            position = _auto_mass_label_position(
                ax,
                curve,
                curves,
                placed_boxes,
                text_size=text_size,
                legend_box_axes=legend_box_axes,
            )

        label_width, label_height = _mass_label_size_data(
            ax,
            str(curve.series.count),
            text_size,
        )
        placed_boxes.append(
            (
                position[0],
                position[1],
                position[0] + label_width,
                position[1] + label_height,
            )
        )
        positions[curve.series.key] = position

    return positions


def plot_mass_by_wavelength(
    data: MassByWavelengthData,
    *,
    ax=None,
    bandwidth: float = 0.5,
    figure_size: tuple[float, float] = FIGURE_SIZE,
    text_size: int = FIGURE_TEXT_SIZE,
    legend_step: float = 0.05,
    label_mode: str = "fixed",
    label_overrides: dict[str, tuple[float, float]] | None = None,
):
    """Plot KDE molecular-mass distributions by detection wavelength."""
    plt = _pyplot()

    created_axes = ax is None
    if ax is None:
        figure = plt.figure(num="Detections at Wavelengths by Mass", figsize=figure_size)
        ax = figure.add_subplot(111)
    else:
        figure = ax.figure

    x_values = np.arange(0, 160)
    curves = []
    for series in data.series:
        if series.count < 2:
            continue
        plot_x = x_values
        if series.truncate_at_data_max:
            max_mass = max(mass for mass in series.masses if mass < x_values[-1])
            plot_x = x_values[: int(max_mass) + 1]
        density = _kde_values(
            series.masses,
            plot_x,
            bandwidth=bandwidth,
        )
        curves.append(
            _MassByWavelengthCurve(
                series=series,
                x_values=plot_x,
                density=density,
            )
        )
        zorder = series.zorder
        ax.plot(
            plot_x,
            density,
            color=series.color,
            zorder=zorder,
        )
        ax.fill_between(
            plot_x,
            density,
            0,
            facecolor=series.color,
            alpha=0.25,
            zorder=zorder,
        )

    header = ax.annotate(
        "Detection Wavelengths",
        xy=(0.95, 0.97),
        xycoords="axes fraction",
        color="black",
        ha="right",
        va="top",
        size=text_size,
    )
    for row, series in enumerate(data.series):
        ax.annotate(
            series.label,
            xy=(0.95, 0.90 - legend_step * row),
            xycoords="axes fraction",
            color=series.color,
            ha="right",
            va="top",
            size=text_size,
        )

    ax.set_xlabel("Atomic Mass (amu)", fontsize=text_size)
    ax.set_ylabel("Probability Density Estimate", fontsize=text_size)
    ax.tick_params(
        axis="both",
        which="both",
        direction="in",
        length=5,
        width=FIGURE_TICK_WIDTH,
        labelsize=text_size,
        colors="black",
    )
    ax.grid(False)
    ax.set_title("")
    if created_axes:
        figure.tight_layout()

    label_positions = _mass_label_positions(
        ax,
        tuple(curves),
        text_size=text_size,
        label_mode=label_mode,
        label_overrides=label_overrides,
        legend_box_axes=(0.52, 0.60, 0.99, 1.0),
    )
    for series in data.series:
        if series.count < 2:
            continue
        ax.annotate(
            str(series.count),
            xy=label_positions[series.key],
            xycoords="data",
            ha="left",
            va="bottom",
            color=series.color,
            size=text_size,
        )
    _underline_axes_text(ax, header)
    return figure, ax


def write_mass_by_wavelength_plot(
    data: MassByWavelengthData,
    output_path: str | Path,
    *,
    bandwidth: float = 0.5,
    label_mode: str = "fixed",
    label_overrides: dict[str, tuple[float, float]] | None = None,
    file_format: str | None = None,
) -> Path:
    """Write the mass-by-wavelength KDE plot and return the output path."""
    output_path = Path(output_path)
    if file_format is None:
        file_format = output_path.suffix.lstrip(".") or "pdf"

    figure, _ = plot_mass_by_wavelength(
        data,
        bandwidth=bandwidth,
        label_mode=label_mode,
        label_overrides=label_overrides,
    )
    figure.savefig(
        output_path,
        format=file_format,
        transparent=True,
        bbox_inches="tight",
    )
    _pyplot().close(figure)
    return output_path


def _stable_jitter_seed(value: str) -> int:
    """Return a deterministic seed for figure jitter."""
    return sum((index + 1) * ord(char) for index, char in enumerate(value))


def _mass_boxplot_color(series: MassByWavelengthSeries) -> str:
    """Return the production color for a mass-by-wavelength boxplot series."""
    return MASS_BY_WAVELENGTH_BOX_COLORS.get(series.key, series.color)


def _source_boxplot_color(category) -> str:
    """Return the production color for a source-type boxplot category."""
    return DU_BY_SOURCE_TYPE_BOX_COLORS.get(category.key, category.color)


def _du_source_boxplot_color(category: DUBySourceTypeCategory) -> str:
    """Return the production color for a DU-by-source boxplot category."""
    return _source_boxplot_color(category)


def plot_mass_by_wavelength_boxplot(
    data: MassByWavelengthData,
    *,
    ax=None,
    figure_size: tuple[float, float] = MASS_BY_WAVELENGTH_BOX_SIZE,
    text_size: float = 18.0,
    tick_size: float = 16.0,
    point_size: float = 22.0,
    point_alpha: float = 0.34,
    x_limit: tuple[float, float] = (0, 240),
    marker_mass: float = 80.0,
):
    """Plot wavelength mass distributions as horizontal boxes plus detections."""
    plt = _pyplot()

    if ax is None:
        figure = plt.figure(
            num="Mass Distributions by Detection Wavelength",
            figsize=figure_size,
        )
        ax = figure.add_subplot(111)
    else:
        figure = ax.figure

    series = tuple(item for item in data.series if item.count > 0)
    y_positions = np.arange(len(series), 0, -1)
    mass_lists = [np.array(item.masses, dtype=float) for item in series]

    boxplot = ax.boxplot(
        mass_lists,
        positions=y_positions,
        vert=False,
        widths=0.42,
        patch_artist=True,
        showfliers=False,
        whis=(10, 90),
        manage_ticks=False,
    )

    for index, item in enumerate(series):
        color = _mass_boxplot_color(item)
        box = boxplot["boxes"][index]
        box.set_facecolor(color)
        box.set_alpha(0.14 if item.key == "IR" else 0.18)
        box.set_edgecolor(color)
        box.set_linewidth(2.0)
        box.set_zorder(3)

        for whisker in boxplot["whiskers"][2 * index : 2 * index + 2]:
            whisker.set_color(color)
            whisker.set_linewidth(1.6)
            whisker.set_alpha(0.9)
            whisker.set_zorder(3)
        for cap in boxplot["caps"][2 * index : 2 * index + 2]:
            cap.set_color(color)
            cap.set_linewidth(1.6)
            cap.set_alpha(0.9)
            cap.set_zorder(3)

    for median in boxplot["medians"]:
        median.set_color("black")
        median.set_linewidth(2.2)
        median.set_zorder(5)

    for y_position, item, masses in zip(y_positions, series, mass_lists):
        rng = np.random.default_rng(_stable_jitter_seed(item.key))
        jitter = rng.uniform(-0.135, 0.135, len(masses))
        ax.scatter(
            masses,
            np.full(len(masses), y_position) + jitter,
            s=point_size,
            color=_mass_boxplot_color(item),
            alpha=point_alpha,
            edgecolors="none",
            zorder=4,
            rasterized=True,
        )

    ax.axvline(
        marker_mass,
        color=MIT_GRAY,
        linestyle="--",
        linewidth=1.3,
        alpha=0.65,
        zorder=1,
    )
    ax.text(
        marker_mass + 1.6,
        y_positions[0] + 0.48,
        f"{int(marker_mass)} amu",
        ha="left",
        va="top",
        fontsize=tick_size,
        color="#555555",
    )

    ax.set_yticks(y_positions)
    ax.set_yticklabels(
        [item.key for item in series],
        fontsize=tick_size,
    )
    ax.set_xlabel("Atomic Mass (amu)", fontsize=text_size)
    ax.set_xlim(*x_limit)
    ax.set_ylim(0.45, len(series) + 0.55)
    ax.tick_params(
        axis="both",
        which="both",
        direction="in",
        length=5,
        width=FIGURE_TICK_WIDTH,
        labelsize=tick_size,
        colors="black",
    )
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
    ax.grid(False)
    ax.set_title("")
    figure.tight_layout()
    return figure, ax


def write_mass_by_wavelength_boxplot(
    data: MassByWavelengthData,
    output_path: str | Path,
    *,
    file_format: str | None = None,
) -> Path:
    """Write the production mass-by-wavelength boxplot and return the path."""
    output_path = Path(output_path)
    if file_format is None:
        file_format = output_path.suffix.lstrip(".") or "pdf"

    figure, _ = plot_mass_by_wavelength_boxplot(data)
    figure.savefig(
        output_path,
        format=file_format,
        transparent=True,
        bbox_inches="tight",
    )
    _pyplot().close(figure)
    return output_path


def plot_du_by_source_type_boxplot(
    data: DUBySourceTypeData,
    *,
    ax=None,
    figure_size: tuple[float, float] = MASS_BY_WAVELENGTH_BOX_SIZE,
    text_size: float = 18.0,
    tick_size: float = 16.0,
    point_size: float = 22.0,
    point_alpha: float = 0.34,
    include_negative_du: bool = False,
    x_limit: tuple[float, float] | None = None,
    source_label_overrides: dict[str, str] | None = None,
):
    """Plot DU distributions by source type as boxes plus exact detections."""
    plt = _pyplot()
    source_label_overrides = source_label_overrides or {}

    if ax is None:
        figure = plt.figure(
            num="DU Distributions by Source Type",
            figsize=figure_size,
        )
        ax = figure.add_subplot(111)
    else:
        figure = ax.figure

    categories_by_key = _category_by_key(data.categories)
    categories = tuple(
        categories_by_key[key]
        for key in DU_BY_SOURCE_TYPE_BOXPLOT_ORDER
        if key in categories_by_key
    )
    value_lists = []
    for category in categories:
        values = [
            value
            for value in category.values
            if include_negative_du or value >= 0
        ]
        value_lists.append(np.array(values, dtype=float))

    categories_with_values = tuple(
        category
        for category, values in zip(categories, value_lists)
        if len(values) > 0
    )
    value_lists = [
        values
        for values in value_lists
        if len(values) > 0
    ]
    if not value_lists:
        raise ValueError("No DU values are available for the requested plot.")

    y_positions = np.arange(len(categories_with_values), 0, -1)

    boxplot = ax.boxplot(
        value_lists,
        positions=y_positions,
        vert=False,
        widths=0.42,
        patch_artist=True,
        showfliers=False,
        whis=(10, 90),
        manage_ticks=False,
    )

    for index, category in enumerate(categories_with_values):
        color = _du_source_boxplot_color(category)
        box = boxplot["boxes"][index]
        box.set_facecolor(color)
        box.set_alpha(0.18)
        box.set_edgecolor(color)
        box.set_linewidth(2.0)
        box.set_zorder(3)

        for whisker in boxplot["whiskers"][2 * index : 2 * index + 2]:
            whisker.set_color(color)
            whisker.set_linewidth(1.6)
            whisker.set_alpha(0.9)
            whisker.set_zorder(3)
        for cap in boxplot["caps"][2 * index : 2 * index + 2]:
            cap.set_color(color)
            cap.set_linewidth(1.6)
            cap.set_alpha(0.9)
            cap.set_zorder(3)

    for median in boxplot["medians"]:
        median.set_color("black")
        median.set_linewidth(2.2)
        median.set_zorder(5)

    for y_position, category, values in zip(
        y_positions,
        categories_with_values,
        value_lists,
    ):
        rng = np.random.default_rng(_stable_jitter_seed(category.key))
        jitter = rng.uniform(-0.135, 0.135, len(values))
        ax.scatter(
            values,
            np.full(len(values), y_position) + jitter,
            s=point_size,
            color=_du_source_boxplot_color(category),
            alpha=point_alpha,
            edgecolors="none",
            zorder=4,
            rasterized=True,
        )

    if x_limit is None:
        min_du = min(float(np.min(values)) for values in value_lists)
        max_du = max(float(np.max(values)) for values in value_lists)
        left_limit = min(-0.2, min_du - 0.2) if include_negative_du else -0.2
        x_limit = (left_limit, max(15.5, max_du + 1.5))
    from matplotlib.transforms import blended_transform_factory

    count_label_transform = blended_transform_factory(ax.transAxes, ax.transData)
    for y_position, values in zip(y_positions, value_lists):
        ax.text(
            1.015,
            y_position,
            f"n={len(values)}",
            ha="left",
            va="center",
            fontsize=tick_size,
            color="#444444",
            transform=count_label_transform,
            clip_on=False,
        )

    ax.set_yticks(y_positions)
    ax.set_yticklabels(
        [
            source_label_overrides.get(category.key, category.label)
            for category in categories_with_values
        ],
        fontsize=tick_size,
    )
    ax.set_xlabel("Degree of Unsaturation", fontsize=text_size)
    ax.set_xlim(*x_limit)
    ax.set_ylim(0.45, len(categories_with_values) + 0.55)
    ax.set_xticks(np.arange(0, int(x_limit[1]) + 1, 2))
    ax.tick_params(
        axis="both",
        which="both",
        direction="in",
        length=5,
        width=FIGURE_TICK_WIDTH,
        labelsize=tick_size,
        colors="black",
    )
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
    ax.grid(False)
    ax.set_title("")
    figure.tight_layout()
    return figure, ax


def write_du_by_source_type_boxplot(
    data: DUBySourceTypeData,
    output_path: str | Path,
    *,
    include_negative_du: bool = False,
    source_label_overrides: dict[str, str] | None = None,
    file_format: str | None = None,
) -> Path:
    """Write the production DU-by-source box/strip plot and return the path."""
    output_path = Path(output_path)
    if file_format is None:
        file_format = output_path.suffix.lstrip(".") or "pdf"

    figure, _ = plot_du_by_source_type_boxplot(
        data,
        include_negative_du=include_negative_du,
        source_label_overrides=source_label_overrides,
    )
    figure.savefig(
        output_path,
        format=file_format,
        transparent=True,
        bbox_inches="tight",
    )
    _pyplot().close(figure)
    return output_path


def plot_relative_du_by_source_type_boxplot(
    data: RelativeDUBySourceTypeData,
    *,
    ax=None,
    figure_size: tuple[float, float] = MASS_BY_WAVELENGTH_BOX_SIZE,
    text_size: float = 18.0,
    tick_size: float = 16.0,
    point_size: float = 22.0,
    point_alpha: float = 0.34,
    include_negative_du: bool = False,
    x_limit: tuple[float, float] | None = None,
    source_label_overrides: dict[str, str] | None = None,
):
    """Plot relative-DU source distributions as boxes plus exact values."""
    plt = _pyplot()
    source_label_overrides = source_label_overrides or {}

    if ax is None:
        figure = plt.figure(
            num="Relative DU Distributions by Source Type",
            figsize=figure_size,
        )
        ax = figure.add_subplot(111)
    else:
        figure = ax.figure

    categories_by_key = _category_by_key(data.categories)
    categories = tuple(
        categories_by_key[key]
        for key in DU_BY_SOURCE_TYPE_BOXPLOT_ORDER
        if key in categories_by_key
    )
    value_lists = []
    for category in categories:
        values = [
            value
            for value in category.values
            if include_negative_du or value >= 0
        ]
        value_lists.append(np.array(values, dtype=float))

    categories_with_values = tuple(
        category
        for category, values in zip(categories, value_lists)
        if len(values) > 0
    )
    value_lists = [
        values
        for values in value_lists
        if len(values) > 0
    ]
    if not value_lists:
        raise ValueError("No relative-DU values are available for the plot.")

    y_positions = np.arange(len(categories_with_values), 0, -1)
    boxplot = ax.boxplot(
        value_lists,
        positions=y_positions,
        vert=False,
        widths=0.42,
        patch_artist=True,
        showfliers=False,
        whis=(10, 90),
        manage_ticks=False,
    )

    for index, category in enumerate(categories_with_values):
        color = _du_source_boxplot_color(category)
        box = boxplot["boxes"][index]
        box.set_facecolor(color)
        box.set_alpha(0.18)
        box.set_edgecolor(color)
        box.set_linewidth(2.0)
        box.set_zorder(3)

        for whisker in boxplot["whiskers"][2 * index : 2 * index + 2]:
            whisker.set_color(color)
            whisker.set_linewidth(1.6)
            whisker.set_alpha(0.9)
            whisker.set_zorder(3)
        for cap in boxplot["caps"][2 * index : 2 * index + 2]:
            cap.set_color(color)
            cap.set_linewidth(1.6)
            cap.set_alpha(0.9)
            cap.set_zorder(3)

    for median in boxplot["medians"]:
        median.set_color("black")
        median.set_linewidth(2.2)
        median.set_zorder(5)

    for y_position, category, values in zip(
        y_positions,
        categories_with_values,
        value_lists,
    ):
        rng = np.random.default_rng(_stable_jitter_seed(category.key))
        jitter = rng.uniform(-0.135, 0.135, len(values))
        ax.scatter(
            values,
            np.full(len(values), y_position) + jitter,
            s=point_size,
            color=_du_source_boxplot_color(category),
            alpha=point_alpha,
            edgecolors="none",
            zorder=4,
            rasterized=True,
        )

    if x_limit is None:
        min_relative_du = min(float(np.min(values)) for values in value_lists)
        left_limit = (
            min(-0.04, min_relative_du - 0.04)
            if include_negative_du
            else -0.04
        )
        x_limit = (left_limit, 1.05)

    from matplotlib.transforms import blended_transform_factory

    count_label_transform = blended_transform_factory(ax.transAxes, ax.transData)
    for y_position, values in zip(y_positions, value_lists):
        ax.text(
            1.015,
            y_position,
            f"n={len(values)}",
            ha="left",
            va="center",
            fontsize=tick_size,
            color="#444444",
            transform=count_label_transform,
            clip_on=False,
        )

    ax.set_yticks(y_positions)
    ax.set_yticklabels(
        [
            source_label_overrides.get(category.key, category.label)
            for category in categories_with_values
        ],
        fontsize=tick_size,
    )
    ax.set_xlabel("Relative Degree of Unsaturation", fontsize=text_size)
    ax.set_xlim(*x_limit)
    ax.set_ylim(0.45, len(categories_with_values) + 0.55)
    ax.set_xticks(np.linspace(0, 1.0, 6))
    ax.tick_params(
        axis="both",
        which="both",
        direction="in",
        length=5,
        width=FIGURE_TICK_WIDTH,
        labelsize=tick_size,
        colors="black",
    )
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
    ax.grid(False)
    ax.set_title("")
    figure.tight_layout()
    return figure, ax


def write_relative_du_by_source_type_boxplot(
    data: RelativeDUBySourceTypeData,
    output_path: str | Path,
    *,
    include_negative_du: bool = False,
    source_label_overrides: dict[str, str] | None = None,
    file_format: str | None = None,
) -> Path:
    """Write the production relative-DU source box/strip plot."""
    output_path = Path(output_path)
    if file_format is None:
        file_format = output_path.suffix.lstrip(".") or "pdf"

    figure, _ = plot_relative_du_by_source_type_boxplot(
        data,
        include_negative_du=include_negative_du,
        source_label_overrides=source_label_overrides,
    )
    figure.savefig(
        output_path,
        format=file_format,
        transparent=True,
        bbox_inches="tight",
    )
    _pyplot().close(figure)
    return output_path


def plot_mass_by_source_type_boxplot(
    data: MassBySourceTypeData,
    *,
    ax=None,
    figure_size: tuple[float, float] = MASS_BY_WAVELENGTH_BOX_SIZE,
    text_size: float = 18.0,
    tick_size: float = 16.0,
    point_size: float = 22.0,
    point_alpha: float = 0.34,
    x_limit: tuple[float, float] | None = None,
    source_label_overrides: dict[str, str] | None = None,
):
    """Plot source-type mass distributions as boxes plus exact detections."""
    plt = _pyplot()
    source_label_overrides = source_label_overrides or {}

    if ax is None:
        figure = plt.figure(
            num="Mass Distributions by Source Type",
            figsize=figure_size,
        )
        ax = figure.add_subplot(111)
    else:
        figure = ax.figure

    categories_by_key = _category_by_key(data.categories)
    categories = tuple(
        categories_by_key[key]
        for key in DU_BY_SOURCE_TYPE_BOXPLOT_ORDER
        if key in categories_by_key
    )
    mass_lists = [
        np.array(category.masses, dtype=float)
        for category in categories
    ]
    categories_with_values = tuple(
        category
        for category, masses in zip(categories, mass_lists)
        if len(masses) > 0
    )
    mass_lists = [
        masses
        for masses in mass_lists
        if len(masses) > 0
    ]
    if not mass_lists:
        raise ValueError("No molecular masses are available for the requested plot.")

    y_positions = np.arange(len(categories_with_values), 0, -1)
    boxplot = ax.boxplot(
        mass_lists,
        positions=y_positions,
        vert=False,
        widths=0.42,
        patch_artist=True,
        showfliers=False,
        whis=(10, 90),
        manage_ticks=False,
    )

    for index, category in enumerate(categories_with_values):
        color = _source_boxplot_color(category)
        box = boxplot["boxes"][index]
        box.set_facecolor(color)
        box.set_alpha(0.18)
        box.set_edgecolor(color)
        box.set_linewidth(2.0)
        box.set_zorder(3)

        for whisker in boxplot["whiskers"][2 * index : 2 * index + 2]:
            whisker.set_color(color)
            whisker.set_linewidth(1.6)
            whisker.set_alpha(0.9)
            whisker.set_zorder(3)
        for cap in boxplot["caps"][2 * index : 2 * index + 2]:
            cap.set_color(color)
            cap.set_linewidth(1.6)
            cap.set_alpha(0.9)
            cap.set_zorder(3)

    for median in boxplot["medians"]:
        median.set_color("black")
        median.set_linewidth(2.2)
        median.set_zorder(5)

    for y_position, category, masses in zip(
        y_positions,
        categories_with_values,
        mass_lists,
    ):
        rng = np.random.default_rng(_stable_jitter_seed(category.key))
        jitter = rng.uniform(-0.135, 0.135, len(masses))
        ax.scatter(
            masses,
            np.full(len(masses), y_position) + jitter,
            s=point_size,
            color=_source_boxplot_color(category),
            alpha=point_alpha,
            edgecolors="none",
            zorder=4,
            rasterized=True,
        )

    if x_limit is None:
        max_mass = max(float(np.max(masses)) for masses in mass_lists)
        x_limit = (0, max(160, int(np.ceil((max_mass + 10) / 20) * 20)))

    from matplotlib.transforms import blended_transform_factory

    count_label_transform = blended_transform_factory(ax.transAxes, ax.transData)
    for y_position, masses in zip(y_positions, mass_lists):
        ax.text(
            1.015,
            y_position,
            f"n={len(masses)}",
            ha="left",
            va="center",
            fontsize=tick_size,
            color="#444444",
            transform=count_label_transform,
            clip_on=False,
        )

    ax.set_yticks(y_positions)
    ax.set_yticklabels(
        [
            source_label_overrides.get(category.key, category.label)
            for category in categories_with_values
        ],
        fontsize=tick_size,
    )
    ax.set_xlabel("Molecular Mass (amu)", fontsize=text_size)
    ax.set_xlim(*x_limit)
    ax.set_ylim(0.45, len(categories_with_values) + 0.55)
    ax.tick_params(
        axis="both",
        which="both",
        direction="in",
        length=5,
        width=FIGURE_TICK_WIDTH,
        labelsize=tick_size,
        colors="black",
    )
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
    ax.grid(False)
    ax.set_title("")
    figure.tight_layout()
    return figure, ax


def write_mass_by_source_type_boxplot(
    data: MassBySourceTypeData,
    output_path: str | Path,
    *,
    source_label_overrides: dict[str, str] | None = None,
    file_format: str | None = None,
) -> Path:
    """Write the production mass-by-source box/strip plot."""
    output_path = Path(output_path)
    if file_format is None:
        file_format = output_path.suffix.lstrip(".") or "pdf"

    figure, _ = plot_mass_by_source_type_boxplot(
        data,
        source_label_overrides=source_label_overrides,
    )
    figure.savefig(
        output_path,
        format=file_format,
        transparent=True,
        bbox_inches="tight",
    )
    _pyplot().close(figure)
    return output_path


def plot_molecules_by_wavelength_atoms(
    data: MoleculesByWavelengthAtomsData,
    *,
    axes=None,
    bandwidth: float = 0.5,
    figure_size: tuple[float, float] = FIGURE_SIZE,
    trace_color: str = ASTROMOL_BLUE,
    text_size: int = 18,
    panel_label_size: int = FIGURE_TEXT_SIZE,
):
    """Plot atom-count distributions grouped by first-detection wavelength."""
    plt = _pyplot()

    if axes is None:
        figure, axes_array = plt.subplots(
            2,
            3,
            num="Molecules Detected in Each Wavelength by Number of Atoms",
            figsize=figure_size,
        )
        axes_list = tuple(axes_array.flat)
    else:
        axes_list = tuple(axes)
        figure = axes_list[0].figure

    x_values = np.arange(0, data.max_atoms + 1, 0.5)
    for index, (axis, series) in enumerate(zip(axes_list, data.series)):
        axis.tick_params(
            axis="both",
            which="both",
            direction="in",
            length=5,
            width=FIGURE_TICK_WIDTH,
            labelsize=text_size,
            colors="black",
        )
        for spine in axis.spines.values():
            spine.set_visible(True)
            spine.set_color("black")

        values = tuple(float(value) for value in series.atom_counts)
        if series.plot == "kde" and len(set(values)) > 1:
            density = _kde_values(values, x_values, bandwidth=bandwidth)
            axis.plot(x_values, density, color=trace_color)
            axis.fill_between(
                x_values,
                density,
                0,
                facecolor=trace_color,
                alpha=0.25,
            )
            axis.set_ylim(0, 0.65)
        else:
            axis.hist(
                series.atom_counts,
                bins=[0.5, 1.5, 2.5, 3.5, 4.5],
                color=trace_color,
                alpha=0.75,
            )
            axis.set_xlim(0, 20)
            axis.set_ylim(0, 8)

        axis.annotate(
            series.label,
            xy=(0.95, 0.96),
            xycoords="axes fraction",
            ha="right",
            va="top",
            size=panel_label_size,
            color="black",
        )
        axis.set_xticks([0, 5, 10, 15, 20])
        if index in {0, 1, 2, 4, 5}:
            axis.set_xticklabels([])

    axes_list[0].set_ylabel("Probability Density Estimate", fontsize=text_size)
    axes_list[3].set_xlabel("# of Atoms", fontsize=text_size)
    axes_list[3].set_ylabel("Probability Density Estimate", fontsize=text_size)
    axes_list[4].set_ylabel("# of Detected Molecules", fontsize=text_size)
    figure.tight_layout()
    return figure, axes_list


def write_molecules_by_wavelength_atoms_plot(
    data: MoleculesByWavelengthAtomsData,
    output_path: str | Path,
    *,
    bandwidth: float = 0.5,
    file_format: str | None = None,
) -> Path:
    """Write the atom-count-by-wavelength plot and return the output path."""
    output_path = Path(output_path)
    if file_format is None:
        file_format = output_path.suffix.lstrip(".") or "pdf"

    figure, _ = plot_molecules_by_wavelength_atoms(
        data,
        bandwidth=bandwidth,
    )
    figure.savefig(
        output_path,
        format=file_format,
        transparent=True,
        bbox_inches="tight",
        pad_inches=0,
    )
    _pyplot().close(figure)
    return output_path


def plot_molecules_by_wavelength_atoms_bubble_heatmap(
    data: MoleculesByWavelengthAtomsData,
    *,
    ax=None,
    colorbar_ax=None,
    value_mode: str = "count",
    figure_size: tuple[float, float] = MOLECULES_BY_WAVELENGTH_BUBBLE_SIZE,
    text_size: float = 18.0,
    tick_size: float = 16.0,
    value_text_size: float = 12.0,
    colorbar_text_size: float = 12.0,
    min_bubble_size: float = 85.0,
    max_bubble_size: float = 1300.0,
    value_text_y_offset: float = 0.025,
):
    """Plot atom-count/wavelength counts as a discrete bubble heatmap."""
    plt = _pyplot()
    if value_mode not in {"count", "row_percent"}:
        msg = "value_mode must be 'count' or 'row_percent'"
        raise ValueError(msg)

    if ax is None:
        figure = plt.figure(
            num="Molecules by Wavelength and Atom Count",
            figsize=figure_size,
        )
        ax = figure.add_axes((0.13, 0.16, 0.76, 0.78))
    else:
        figure = ax.figure

    if colorbar_ax is None:
        colorbar_ax = figure.add_axes((0.92, 0.16, 0.035, 0.78))

    bins = MOLECULES_BY_WAVELENGTH_ATOM_BINS
    bin_labels = [str(item) for item in bins]
    row_labels = [series.key for series in data.series]
    count_matrix = data.matrix(bins)
    if value_mode == "row_percent":
        row_totals = count_matrix.sum(axis=1, keepdims=True)
        plot_matrix = np.divide(
            count_matrix,
            row_totals,
            out=np.zeros_like(count_matrix, dtype=float),
            where=row_totals != 0,
        ) * 100
        colorbar_label = "Row Percentage (%)"
    else:
        plot_matrix = count_matrix.astype(float)
        colorbar_label = "Molecules"

    max_value = max(float(plot_matrix.max()), 1.0)
    effective_min_bubble_size = min_bubble_size
    if value_mode == "row_percent":
        effective_min_bubble_size = max(effective_min_bubble_size, 165.0)
    colormap = _periodic_heatmap_colormap()

    from matplotlib import colors

    norm = colors.Normalize(vmin=0, vmax=max_value)

    for row in range(plot_matrix.shape[0]):
        for col in range(plot_matrix.shape[1]):
            count = int(count_matrix[row, col])
            plot_value = float(plot_matrix[row, col])
            if count == 0:
                continue
            bubble_size = (
                effective_min_bubble_size
                + (max_bubble_size - effective_min_bubble_size)
                * (plot_value / max_value) ** 0.72
            )
            facecolor = colormap(norm(plot_value))
            text_value = f"{plot_value:.0f}" if value_mode == "row_percent" else str(count)
            ax.scatter(
                col,
                row,
                s=bubble_size,
                facecolor=facecolor,
                edgecolor="black",
                linewidth=0.7,
                alpha=0.88,
                zorder=3,
            )
            ax.text(
                col,
                row + value_text_y_offset,
                text_value,
                ha="center",
                va="center",
                fontsize=value_text_size,
                color="black",
                zorder=4,
            )

    ax.set_xlim(-0.6, len(bins) - 0.4)
    ax.set_ylim(len(row_labels) - 0.5, -0.5)
    ax.set_xticks(np.arange(len(bins)))
    ax.set_xticklabels(bin_labels, fontsize=tick_size)
    ax.set_yticks(np.arange(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=tick_size)
    ax.set_xlabel("# of Atoms", fontsize=text_size)
    ax.tick_params(
        axis="both",
        which="both",
        direction="in",
        length=5,
        width=FIGURE_TICK_WIDTH,
        colors="black",
    )
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
    ax.grid(False)
    ax.set_title("")

    scalar_mappable = plt.cm.ScalarMappable(cmap=colormap, norm=norm)
    scalar_mappable.set_array([])
    colorbar = figure.colorbar(scalar_mappable, cax=colorbar_ax)
    colorbar.set_label(colorbar_label, fontsize=colorbar_text_size)
    colorbar.ax.tick_params(labelsize=colorbar_text_size, colors="black")

    return figure, ax


def write_molecules_by_wavelength_atoms_bubble_heatmap(
    data: MoleculesByWavelengthAtomsData,
    output_path: str | Path,
    *,
    file_format: str | None = None,
    value_mode: str = "count",
) -> Path:
    """Write the production atom-count/wavelength bubble heatmap."""
    output_path = Path(output_path)
    if file_format is None:
        file_format = output_path.suffix.lstrip(".") or "pdf"

    figure, _ = plot_molecules_by_wavelength_atoms_bubble_heatmap(
        data,
        value_mode=value_mode,
    )
    figure.savefig(
        output_path,
        format=file_format,
        transparent=True,
        bbox_inches="tight",
    )
    _pyplot().close(figure)
    return output_path


def plot_rolling_rate_by_atoms_heatmap(
    data: RollingRateHeatmapData,
    *,
    ax=None,
    colorbar_ax=None,
    cmap: str = "cividis",
    text_size: float = 8.5,
    tick_size: float = 7.2,
    colorbar_text_size: float = 7.2,
):
    """Plot rolling detection rates by atom-count category as a heatmap."""
    plt = _pyplot()
    if ax is None:
        figure = plt.figure(
            num="Rolling Detection Rate By Atoms",
            figsize=ROLLING_RATE_HEATMAP_SIZE,
        )
        ax = figure.add_axes(ROLLING_RATE_HEATMAP_AXES_BOUNDS)
    else:
        figure = ax.figure

    if colorbar_ax is None:
        colorbar_ax = figure.add_axes(ROLLING_RATE_HEATMAP_COLORBAR_BOUNDS)

    image = ax.imshow(
        data.matrix,
        aspect="auto",
        cmap=cmap,
        extent=(
            data.start_year - 0.5,
            data.end_year + 0.5,
            len(data.labels) - 0.5,
            -0.5,
        ),
        interpolation="nearest",
        origin="upper",
        vmin=0,
        vmax=data.vmax,
    )

    ax.set_xlabel("Year", fontsize=text_size)
    ax.set_ylabel("# Atoms", fontsize=text_size, labelpad=1)
    ax.set_yticks(np.arange(len(data.labels)))
    ax.set_yticklabels(data.labels, fontsize=tick_size)
    ax.set_xticks(
        [
            tick
            for tick in ROLLING_RATE_HEATMAP_XTICKS
            if data.start_year <= tick <= data.end_year
        ]
    )
    ax.tick_params(
        axis="x",
        direction="out",
        length=3,
        width=0.8,
        labelsize=tick_size,
        pad=1.5,
    )
    ax.tick_params(
        axis="y",
        direction="out",
        length=0,
        width=0.8,
        labelsize=tick_size,
        pad=2,
    )

    colorbar = figure.colorbar(
        image,
        cax=colorbar_ax,
        orientation="horizontal",
    )
    colorbar.ax.xaxis.set_ticks_position("top")
    colorbar.ax.xaxis.set_label_position("top")
    colorbar.ax.tick_params(labelsize=colorbar_text_size, pad=1.5)
    colorbar.set_label(
        f"{data.window}-yr rolling detections yr$^{{-1}}$",
        fontsize=colorbar_text_size,
        labelpad=2,
    )

    return figure, ax


def write_rolling_rate_by_atoms_heatmap(
    data: RollingRateHeatmapData,
    output_path: str | Path,
    *,
    file_format: str | None = None,
) -> Path:
    """Write the rolling-rate-by-atoms heatmap and return the output path."""
    output_path = Path(output_path)
    if file_format is None:
        file_format = output_path.suffix.lstrip(".") or "pdf"

    figure, _ = plot_rolling_rate_by_atoms_heatmap(data)
    figure.savefig(
        output_path,
        format=file_format,
        transparent=True,
    )
    _pyplot().close(figure)
    return output_path


def _marker_sizes_by_count(
    points: tuple[DetectionRateByAtomsPoint, ...],
    *,
    marker_scale: float = 200,
) -> np.ndarray:
    """Return legacy marker areas proportional to category counts."""
    counts = np.array(
        [
            point.count if point.count > 0 else np.nan
            for point in points
        ],
        dtype=float,
    )
    if np.all(np.isnan(counts)):
        return counts
    return counts * marker_scale / np.nanmin(counts)


def _absolute_marker_sizes_by_count(
    points: tuple[DetectionRateByAtomsPoint, ...],
    *,
    marker_area_per_detection: float = 70,
) -> np.ndarray:
    """Return marker areas on an absolute detections-to-area scale."""
    return np.array(
        [
            point.count * marker_area_per_detection
            if point.count > 0
            else np.nan
            for point in points
        ],
        dtype=float,
    )


def _visible_rate_points(
    data: DetectionRateByAtomsData,
    xlim: tuple[float, float],
) -> tuple[list[int], tuple[DetectionRateByAtomsPoint, ...]]:
    """Return indices and points visible within a rate-by-atoms x-axis range."""
    visible_indices = [
        index
        for index, point in enumerate(data.points)
        if xlim[0] <= point.x_position <= xlim[1]
    ]
    return visible_indices, tuple(data.points[index] for index in visible_indices)


def plot_detection_rate_by_atoms(
    data: DetectionRateByAtomsData,
    *,
    ax=None,
    facecolor: str = ASTROMOL_BLUE,
    edgecolor: str = NRAO_BLUE,
    alpha: float = 0.9,
    text_size: int = FIGURE_TEXT_SIZE,
    xlim: tuple[float, float] = (1, 12.8),
    ylim: tuple[float, float] | None = None,
):
    """Plot average detections/year by atom-count category."""
    plt = _pyplot()
    if ax is None:
        figure, ax = plt.subplots(
            num="Detects Per Year Per Atom",
            figsize=FIGURE_SIZE,
        )
    else:
        figure = ax.figure

    visible_indices, visible_points = _visible_rate_points(data, xlim)
    x_positions = np.array([point.x_position for point in visible_points])
    rates = np.array([point.rate for point in visible_points])
    marker_sizes = _marker_sizes_by_count(data.points)
    visible_marker_sizes = marker_sizes[visible_indices]

    ax.scatter(
        x_positions,
        rates,
        marker="o",
        c=facecolor,
        edgecolors=edgecolor,
        s=visible_marker_sizes,
        alpha=alpha,
    )

    _style_manuscript_axes(
        ax,
        xlabel="Number of Atoms",
        ylabel="Detections/Year*",
        text_size=text_size,
    )
    ax.set_xticks(DETECTION_RATE_BY_ATOMS_XTICKS)
    ax.set_xticklabels(DETECTION_RATE_BY_ATOMS_XTICK_LABELS)
    ax.set_xlim(xlim)
    if ylim is None:
        visible_rates = [
            point.rate
            for point in data.visible_points
            if not np.isnan(point.rate)
        ]
        upper_limit = max(1.0, max(visible_rates) * 1.1)
        ylim = (0, upper_limit)
    ax.set_ylim(ylim)

    ax.annotate(
        "*Since year of first detection",
        xy=(0.35, 0.9),
        xycoords="axes fraction",
        ha="left",
        size=text_size,
    )
    ax.annotate(
        "Marker size proportional to\n"
        "total # of detections",
        xy=(0.05, 0.1),
        xycoords="axes fraction",
        ha="left",
        size=text_size,
    )
    for label in ax.get_xmajorticklabels():
        if label.get_text() in {"Fullerenes", "PAHs"}:
            label.set_rotation(-45)

    _finalize_manuscript_figure(figure, ax)
    return figure, ax


def write_detection_rate_by_atoms_plot(
    data: DetectionRateByAtomsData,
    output_path: str | Path,
    *,
    file_format: str | None = None,
) -> Path:
    """Write the detection-rate-by-atoms plot and return the output path."""
    output_path = Path(output_path)
    if file_format is None:
        file_format = output_path.suffix.lstrip(".") or "pdf"

    figure, _ = plot_detection_rate_by_atoms(data)
    figure.savefig(
        output_path,
        format=file_format,
        transparent=True,
        bbox_inches="tight",
    )
    _pyplot().close(figure)
    return output_path


def plot_detection_rate_by_atoms_comparison(
    current_data: DetectionRateByAtomsData,
    baseline_data: DetectionRateByAtomsData,
    *,
    ax=None,
    current_label: str = "Current",
    baseline_label: str = "2021",
    marker_area_per_detection: float = 70,
    current_facecolor: str = ASTROMOL_BLUE,
    current_edgecolor: str = NRAO_BLUE,
    baseline_edgecolor: str = "black",
    current_alpha: float = 0.86,
    text_size: int = FIGURE_TEXT_SIZE,
    legend_text_size: int = 18,
    xlim: tuple[float, float] = (1, 12.8),
    ylim: tuple[float, float] = (0, 1.0),
):
    """Plot current rate-by-atoms values over a 2021 comparison baseline.

    Marker area uses a single absolute scale for both datasets, so area is
    directly comparable between the baseline and current views.
    """
    plt = _pyplot()
    if ax is None:
        figure, ax = plt.subplots(
            num="Detection Rate By Atoms Comparison",
            figsize=FIGURE_SIZE,
        )
    else:
        figure = ax.figure

    baseline_indices, baseline_points = _visible_rate_points(baseline_data, xlim)
    current_indices, current_points = _visible_rate_points(current_data, xlim)
    baseline_sizes = _absolute_marker_sizes_by_count(
        baseline_data.points,
        marker_area_per_detection=marker_area_per_detection,
    )[baseline_indices]
    current_sizes = _absolute_marker_sizes_by_count(
        current_data.points,
        marker_area_per_detection=marker_area_per_detection,
    )[current_indices]

    baseline = ax.scatter(
        [point.x_position for point in baseline_points],
        [point.rate for point in baseline_points],
        s=baseline_sizes,
        marker="o",
        facecolors="none",
        edgecolors=baseline_edgecolor,
        linewidths=1.0,
        zorder=1,
    )
    baseline.set_hatch("///")

    ax.scatter(
        [point.x_position for point in current_points],
        [point.rate for point in current_points],
        s=current_sizes,
        marker="o",
        c=current_facecolor,
        edgecolors=current_edgecolor,
        linewidths=1.3,
        alpha=current_alpha,
        zorder=2,
    )

    _style_manuscript_axes(
        ax,
        xlabel="Number of Atoms",
        ylabel="Detections/Year*",
        text_size=text_size,
    )
    ax.set_xticks(DETECTION_RATE_BY_ATOMS_XTICKS)
    ax.set_xticklabels(DETECTION_RATE_BY_ATOMS_XTICK_LABELS)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    ax.annotate(
        "*Since year of first detection",
        xy=(0.35, 0.9),
        xycoords="axes fraction",
        ha="left",
        size=text_size,
    )
    ax.annotate(
        "Marker area proportional to\n"
        "total # of detections",
        xy=(0.05, 0.1),
        xycoords="axes fraction",
        ha="left",
        size=text_size,
    )

    from matplotlib.lines import Line2D

    legend_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor=current_facecolor,
            markeredgecolor=current_edgecolor,
            markersize=15,
            label=current_label,
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color=baseline_edgecolor,
            markerfacecolor="none",
            markeredgecolor=baseline_edgecolor,
            markersize=15,
            label=baseline_label,
        ),
    ]
    ax.legend(
        handles=legend_handles,
        loc="upper right",
        bbox_to_anchor=(0.98, 0.78),
        frameon=False,
        fontsize=legend_text_size,
    )

    _finalize_manuscript_figure(figure, ax)
    return figure, ax


def write_detection_rate_by_atoms_comparison_plot(
    current_data: DetectionRateByAtomsData,
    baseline_data: DetectionRateByAtomsData,
    output_path: str | Path,
    *,
    file_format: str | None = None,
    current_label: str = "Current",
    baseline_label: str = "2021",
) -> Path:
    """Write the production detection-rate-by-atoms comparison plot."""
    output_path = Path(output_path)
    if file_format is None:
        file_format = output_path.suffix.lstrip(".") or "pdf"

    figure, _ = plot_detection_rate_by_atoms_comparison(
        current_data,
        baseline_data,
        current_label=current_label,
        baseline_label=baseline_label,
    )
    figure.savefig(
        output_path,
        format=file_format,
        transparent=True,
    )
    _pyplot().close(figure)
    return output_path


def plot_facility_shares(
    data: FacilityShareData,
    *,
    axes=None,
    active_color: str = ASTROMOL_BLUE,
    inactive_color: str = FACILITY_SHARE_INACTIVE,
    remainder_color: str = FACILITY_SHARE_BACKGROUND,
    title_text_size: int = 17,
    date_text_size: int = 14,
    percent_text_size: int = 15,
    figure_size: tuple[float, float] = (8, 8),
    pie_radius: float = 0.9,
):
    """Plot facility shares as a 3x3 pie grid."""
    plt = _pyplot()
    import matplotlib.patches as patches
    import numpy as np

    if axes is None:
        figure, axes = plt.subplots(
            3,
            3,
            num="Facility Shares",
            figsize=figure_size,
        )
    else:
        figure = axes.flat[0].figure

    flat_axes = list(axes.flat)
    for index, ax in enumerate(flat_axes):
        if index >= len(data.facilities):
            ax.axis("off")
            continue

        facility = data.facilities[index]
        facility_color = active_color if facility.active_at_view else inactive_color
        fraction = facility.fraction
        pie_values = [fraction, 1.0 - fraction]

        if index > 5:
            slices, _ = ax.pie(
                pie_values,
                colors=[facility_color, remainder_color],
                wedgeprops={"linewidth": 1.0, "edgecolor": "black"},
                radius=pie_radius,
            )
            angle = (slices[0].theta2 - slices[0].theta1) / 2.0 + slices[0].theta1
            y_position = pie_radius * np.sin(np.deg2rad(angle))
            x_position = pie_radius * np.cos(np.deg2rad(angle))
            horizontal_alignment = {-1: "right", 1: "left"}[
                int(np.sign(x_position)) or 1
            ]
            ax.annotate(
                f"{facility.percent}%",
                xy=(x_position, y_position),
                xytext=(1.01 * np.sign(x_position or 1), y_position),
                horizontalalignment=horizontal_alignment,
                size=percent_text_size,
                fontweight="bold",
                arrowprops={
                    "arrowstyle": "-",
                    "connectionstyle": f"angle,angleA=0,angleB={angle}",
                },
                zorder=5,
                va="center",
            )
        else:
            slices, labels = ax.pie(
                pie_values,
                labels=[f"{facility.percent}%", ""],
                colors=[facility_color, remainder_color],
                wedgeprops={"linewidth": 1.0, "edgecolor": "black"},
                labeldistance=0.7,
                radius=pie_radius,
            )
            for label in labels:
                label.set_horizontalalignment("center")
                label.set_color("white")
                label.set_fontsize(percent_text_size)
                label.set_weight("bold")

        ax.annotate(
            facility.label,
            xy=(0.5, 1.02),
            xycoords="axes fraction",
            ha="center",
            size=title_text_size,
            fontweight="normal",
        )
        ax.annotate(
            f"{facility.start_year} - {facility.end_year}",
            xy=(0.5, 0.925),
            xycoords="axes fraction",
            ha="center",
            size=date_text_size,
            fontweight="normal",
        )

        center = slices[0].center
        radius = slices[0].r
        ax.add_patch(
            patches.Circle(
                center,
                radius,
                fill=False,
                edgecolor="black",
                linewidth=1,
            )
        )

    figure.tight_layout()
    figure.subplots_adjust(wspace=0.0, hspace=0.12)
    return figure, axes


def write_facility_shares_plot(
    data: FacilityShareData,
    output_path: str | Path,
    *,
    file_format: str | None = None,
) -> Path:
    """Write the facility-shares plot and return the output path."""
    output_path = Path(output_path)
    if file_format is None:
        file_format = output_path.suffix.lstrip(".") or "pdf"

    figure, _ = plot_facility_shares(data)
    figure.savefig(
        output_path,
        format=file_format,
        transparent=True,
        bbox_inches="tight",
    )
    _pyplot().close(figure)
    return output_path


def plot_facility_share_bars(
    data: FacilityShareData,
    *,
    ax=None,
    active_color: str = ASTROMOL_BLUE,
    inactive_color: str = FACILITY_SHARE_INACTIVE,
    figure_size: tuple[float, float] = FACILITY_SHARE_BAR_SIZE,
    axes_bounds: tuple[float, float, float, float] = FACILITY_SHARE_BAR_AXES_BOUNDS,
    label_text_size: int = 10,
    value_text_size: int = 10,
    axis_text_size: int = 10,
    tick_text_size: int = 9,
    legend_text_size: int = 9,
):
    """Plot facility shares as a horizontal bar chart.

    This is the production-facing facility-share view. The pie-grid function is
    retained for legacy reproduction, while this bar view is easier to compare
    quantitatively and labels each denominator explicitly.
    """
    plt = _pyplot()

    if ax is None:
        figure = plt.figure(num="Facility Share Bars", figsize=figure_size)
        ax = figure.add_axes(axes_bounds)
    else:
        figure = ax.figure

    facilities = tuple(reversed(data.facilities))
    y_positions = np.arange(len(facilities))
    percents = [facility.percent for facility in facilities]
    colors = [
        active_color if facility.active_at_view else inactive_color
        for facility in facilities
    ]

    ax.barh(
        y_positions,
        percents,
        color=colors,
        edgecolor="black",
        linewidth=0.8,
    )

    for y_position, facility in zip(y_positions, facilities):
        ax.text(
            -0.006,
            y_position + 0.16,
            facility.label,
            transform=ax.get_yaxis_transform(),
            va="center",
            ha="right",
            fontsize=label_text_size,
            fontweight="bold",
            color="black",
            clip_on=False,
        )
        ax.text(
            -0.006,
            y_position - 0.16,
            f"{facility.start_year} - {facility.end_year}",
            transform=ax.get_yaxis_transform(),
            va="center",
            ha="right",
            fontsize=label_text_size,
            fontstyle="italic",
            color="black",
            clip_on=False,
        )
        ax.text(
            facility.percent + 1.0,
            y_position,
            (
                f"{facility.percent}% "
                f"({facility.detection_count}/{facility.total_window_detections})"
            ),
            va="center",
            ha="left",
            fontsize=value_text_size,
            fontweight="bold",
            color="black",
        )

    ax.set_yticks(y_positions)
    ax.set_yticklabels([])
    ax.set_xlabel(
        "Contribution share during facility lifetime (%)",
        fontsize=axis_text_size,
    )
    ax.set_xlim(0, max(74, max(percents) + 12))
    ax.set_ylim(-0.5, len(facilities) - 0.5)
    ax.tick_params(axis="x", labelsize=tick_text_size, length=5, width=1)
    ax.tick_params(axis="y", length=0, pad=2)
    ax.grid(False)
    ax.set_title("")

    plt = _pyplot()
    active_patch = plt.Rectangle(
        (0, 0),
        1,
        1,
        facecolor=active_color,
        edgecolor="black",
        linewidth=0.8,
    )
    inactive_patch = plt.Rectangle(
        (0, 0),
        1,
        1,
        facecolor=inactive_color,
        edgecolor="black",
        linewidth=0.8,
    )
    ax.legend(
        [active_patch, inactive_patch],
        ["Active", "Decommissioned"],
        loc="lower right",
        frameon=False,
        fontsize=legend_text_size,
        handlelength=1.5,
        handleheight=2.4,
        borderpad=0.2,
        labelspacing=0.4,
    )

    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(1.0)
    ax.tick_params(axis="both", colors="black")

    return figure, ax


def write_facility_share_bars_plot(
    data: FacilityShareData,
    output_path: str | Path,
    *,
    file_format: str | None = None,
) -> Path:
    """Write the facility-share bar chart and return the output path."""
    output_path = Path(output_path)
    if file_format is None:
        file_format = output_path.suffix.lstrip(".") or "pdf"

    figure, _ = plot_facility_share_bars(data)
    figure.savefig(
        output_path,
        format=file_format,
        transparent=True,
    )
    _pyplot().close(figure)
    return output_path


def _last_scope_contribution_year(series: ScopeDetectionSeries) -> int:
    """Return the last year where a facility's cumulative count increased."""
    increases = np.where(np.diff(series.counts, prepend=0) > 0)[0]
    if len(increases) == 0:
        return int(series.years[0])
    return int(series.years[increases[-1]])


def _scope_color(series: ScopeDetectionSeries, style: str) -> str:
    """Return the plotting color for a facility series."""
    if style == "modern":
        return MODERN_SCOPE_COLORS_BY_LABEL.get(series.label, series.color)
    return series.color


def _plot_scope_series(
    ax,
    series: ScopeDetectionSeries,
    *,
    style: str,
    line_width: float,
) -> None:
    """Plot one facility contribution series."""
    color = _scope_color(series, style)

    if style == "modern" and series.label == MODERN_SCOPE_HIGHLIGHT_LABEL:
        first_nonzero = np.where(series.counts > 0)[0]
        if len(first_nonzero) == 0:
            return
        start_index = int(first_nonzero[0])
        ax.plot(
            series.years[start_index:],
            series.counts[start_index:],
            color=color,
            linewidth=4.0,
            alpha=1.0,
            solid_capstyle="round",
        )
        return

    if style == "modern" and series.label in MODERN_SCOPE_DORMANT_LABELS:
        last_year = _last_scope_contribution_year(series)
        solid_mask = series.years <= last_year
        dotted_mask = series.years >= last_year
        ax.plot(
            series.years[solid_mask],
            series.counts[solid_mask],
            color=color,
            linewidth=2.15,
            alpha=0.86,
            solid_capstyle="round",
        )
        ax.plot(
            series.years[dotted_mask],
            series.counts[dotted_mask],
            color=color,
            linewidth=2.15,
            alpha=0.50,
            linestyle=(0, (1.2, 2.0)),
            solid_capstyle="round",
        )
        return

    ax.plot(
        series.years,
        series.counts,
        color=color,
        linewidth=2.15 if style == "modern" else line_width,
        alpha=0.86 if style == "modern" else 1.0,
        solid_capstyle="round",
    )


def plot_scopes_by_year(
    data: ScopesByYearData,
    *,
    ax=None,
    figure_size: tuple[float, float] = FIGURE_SIZE,
    axes_bounds: tuple[float, float, float, float] = FIGURE_AXES_BOUNDS,
    text_size: int = FIGURE_TEXT_SIZE,
    annotation_text_size: int = 18,
    line_width: float = 2.0,
    style: str = "legacy",
):
    """Plot cumulative first-detection contributions for prolific facilities."""
    if style not in {"legacy", "modern"}:
        raise ValueError("style must be 'legacy' or 'modern'.")
    plt = _pyplot()

    if ax is None:
        figure = plt.figure(num="Detections Per Facility Over Time", figsize=figure_size)
        ax = figure.add_axes(axes_bounds)
    else:
        figure = ax.figure

    for series in data.series:
        _plot_scope_series(
            ax,
            series,
            style=style,
            line_width=line_width,
        )

    annotation_series = sorted(data.series, key=lambda item: item.rate, reverse=True)
    for row, series in enumerate(annotation_series):
        y_position = 0.90 - 0.05 * row
        color = _scope_color(series, style)
        weight = (
            "bold"
            if style == "modern" and series.label == MODERN_SCOPE_HIGHLIGHT_LABEL
            else "normal"
        )
        ax.annotate(
            series.label,
            xy=(0.10, y_position),
            xycoords="axes fraction",
            color=color,
            ha="left",
            va="top",
            size=annotation_text_size,
            weight=weight,
        )
        ax.annotate(
            f"{series.rate:.1f}/yr",
            xy=(0.37, y_position),
            xycoords="axes fraction",
            color=color,
            ha="left",
            va="top",
            size=annotation_text_size,
            weight=weight,
        )
        ax.annotate(
            f"({series.start_year} - {series.fit_stop_year})",
            xy=(0.47, y_position),
            xycoords="axes fraction",
            color=color,
            ha="left",
            va="top",
            size=annotation_text_size,
            weight=weight,
        )

    ax.set_xlabel("Year", fontsize=text_size)
    ax.set_ylabel("Cumulative Number of Detected Molecules", fontsize=text_size)
    ax.set_xlim(1965, data.end_year + 2)
    ax.set_ylim(0, data.max_count + 5)
    ax.set_xticks([1970, 1980, 1990, 2000, 2010, 2020])
    ax.tick_params(
        axis="both",
        which="both",
        direction="in",
        length=FIGURE_TICK_LENGTH,
        width=FIGURE_TICK_WIDTH,
        labelsize=text_size,
        colors="black",
        top=True,
        right=True,
    )
    ax.yaxis.set_ticks_position("both")
    ax.xaxis.set_ticks_position("both")
    ax.grid(False)
    ax.set_title("")
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(1.0)

    return figure, ax


def write_scopes_by_year_plot(
    data: ScopesByYearData,
    output_path: str | Path,
    *,
    file_format: str | None = None,
    style: str = "legacy",
) -> Path:
    """Write the scopes-by-year plot and return the output path."""
    output_path = Path(output_path)
    if file_format is None:
        file_format = output_path.suffix.lstrip(".") or "pdf"

    figure, _ = plot_scopes_by_year(data, style=style)
    figure.savefig(
        output_path,
        format=file_format,
        transparent=True,
    )
    _pyplot().close(figure)
    return output_path
