"""PowerPoint slide generation helpers for census products.

The slide code is intentionally split into selection, layout, and rendering
steps. That makes the visually tuned PowerPoint output auditable before it is
written to a binary ``.pptx`` file.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import date, datetime
from importlib import metadata
from math import ceil
from pathlib import Path
import re
import subprocess

from .census import CensusView
from .models import Detection, Molecule


DEFAULT_DETECTION_TYPE = "ISM/CSM"
PPD_DETECTION_TYPE = "ppd"
LEGACY_SLIDE_WIDTH_PT = 1920
LEGACY_SLIDE_HEIGHT_PT = 1080
LEGACY_SLIDE_WIDTH_IN = LEGACY_SLIDE_WIDTH_PT / 72
LEGACY_SLIDE_HEIGHT_IN = LEGACY_SLIDE_HEIGHT_PT / 72
LEGACY_MOLECULE_FONT_PT = 22
LEGACY_ROW_SPACING = 1.05
LEGACY_TITLE = "Known Interstellar Molecules"
PPD_TITLE = "Protoplanetary Disk Molecules"
LEGACY_URL = "bmcguir2.github.io/astromol"
LEGACY_CITATION = "McGuire 2022 ApJS 259, 30"
BALANCED_MOLECULE_FONT_CANDIDATES = (22, 21, 20, 19, 18)
COMPACT_MOLECULE_FONT_CANDIDATES = (26, 25, 24, 23, 22)
BALANCED_CONTENT_LEFT = 0.35
BALANCED_CONTENT_RIGHT = 26.62
BALANCED_LABEL_TOP = 1.55
BALANCED_COLUMN_TOP = 2.08
BALANCED_COLUMN_BOTTOM = 14.5
BALANCED_TOP_LABEL_TOP = 1.55
BALANCED_TOP_COLUMN_TOP = 2.08
BALANCED_TOP_COLUMN_BOTTOM = 9.95
BALANCED_BOTTOM_LABEL_TOP = 10.45
BALANCED_BOTTOM_COLUMN_TOP = 10.98
BALANCED_BOTTOM_COLUMN_BOTTOM = 14.5
BALANCED_BOTTOM_LEFT = 15.55
BALANCED_GROUP_GAP = 0.08
BALANCED_COLUMN_GAP = 0.08
BALANCED_MAX_GROUP_GAP = 0.46
BALANCED_FORMULA_WIDTH_FACTOR = 0.010
COMPACT_CONTENT_LEFT = 0.8
COMPACT_CONTENT_RIGHT = 25.9
COMPACT_LABEL_TOP = 2.0
COMPACT_COLUMN_TOP = 2.6
COMPACT_COLUMN_BOTTOM = 11.2
COMPACT_GROUP_GAP = 0.16
COMPACT_MAX_GROUP_GAP = 1.2
BALANCED_TOP_NATOMS = tuple(range(2, 11))
BALANCED_BOTTOM_NATOMS = (11, 12, 13)
BALANCED_TOP_MIN_COLUMNS = {
    2: 2,
    3: 2,
    4: 2,
    5: 2,
    6: 2,
}

FORMULA_TOKEN = re.compile(
    r"(\\ce\{|\\mathrm\{|\\text\{|[{}]|"
    r"\^\{\d+\}|\[\d+[A-Z][a-z]?\]|"
    r"[A-Z][a-z]?|[a-z]+|\d+|[+\-]|.)"
)


@dataclass(frozen=True)
class SlideBox:
    """A PowerPoint box expressed in slide inches."""

    left: float
    top: float
    width: float
    height: float

    @property
    def right(self) -> float:
        """Right edge in inches."""
        return self.left + self.width

    @property
    def bottom(self) -> float:
        """Bottom edge in inches."""
        return self.top + self.height


@dataclass(frozen=True)
class FormulaRun:
    """One formatted PowerPoint text run for a molecule formula."""

    text: str
    baseline: str = "normal"
    italic: bool = False


@dataclass(frozen=True)
class SlideMoleculeEntry:
    """A molecule selected for a slide, with display metadata."""

    molecule: Molecule
    display_formula: str

    @property
    def label(self) -> str:
        """Stable molecule label."""
        return self.molecule.label

    @property
    def natoms(self) -> int:
        """Number of atoms in the molecule."""
        return self.molecule.natoms


@dataclass(frozen=True)
class SlideGroupSpec:
    """Static layout specification for one atom-count group."""

    label: str
    natoms: int
    ncols: int
    label_box: SlideBox
    column_boxes: tuple[SlideBox, ...]
    natoms_greater: bool = False


@dataclass(frozen=True)
class SlideColumnPlan:
    """Molecules assigned to one rendered column."""

    box: SlideBox
    molecules: tuple[SlideMoleculeEntry, ...]
    estimated_capacity: int

    @property
    def overflow(self) -> int:
        """Estimated number of rows beyond the column capacity."""
        return max(0, len(self.molecules) - self.estimated_capacity)


@dataclass(frozen=True)
class SlideGroupPlan:
    """Rendered plan for one atom-count group."""

    spec: SlideGroupSpec
    columns: tuple[SlideColumnPlan, ...]

    @property
    def molecules(self) -> tuple[SlideMoleculeEntry, ...]:
        """All molecules in this group."""
        return tuple(molecule for column in self.columns for molecule in column.molecules)


@dataclass(frozen=True)
class MoleculeSlideLayout:
    """Complete molecule-slide layout plan."""

    groups: tuple[SlideGroupPlan, ...]
    title: str
    total: int
    profile: str = "legacy"
    slide_width_pt: int = LEGACY_SLIDE_WIDTH_PT
    slide_height_pt: int = LEGACY_SLIDE_HEIGHT_PT
    title_box: SlideBox = SlideBox(0.35, 0.22, 18.33, 1.18)
    credit_box: SlideBox = SlideBox(22.15, 0.22, 300 / 72, 90 / 72)
    total_box: SlideBox = SlideBox(6.35, 10.6, 6.0, 1.4)
    footer_box: SlideBox = SlideBox(6.35, 12.2, 6.0, 1.0)
    molecule_font_pt: float = LEGACY_MOLECULE_FONT_PT
    group_label_font_pt: float = 26
    total_font_pt: float = 45
    footer_font_pt: float = 32
    total_line_width_pt: float = 6
    package_version: str = "development"
    last_updated: str | None = None
    detection_type: str = DEFAULT_DETECTION_TYPE
    include_tentative: bool = False
    include_disputed: bool = False
    include_isotopologues: bool = False
    warnings: tuple[str, ...] = field(default_factory=tuple)

    @property
    def molecules(self) -> tuple[SlideMoleculeEntry, ...]:
        """All molecule entries in slide order."""
        return tuple(molecule for group in self.groups for molecule in group.molecules)


def legacy_molecule_slide_specs() -> tuple[SlideGroupSpec, ...]:
    """Return the historically tuned molecule-slide layout specification."""
    return (
        SlideGroupSpec(
            label="2 Atoms",
            natoms=2,
            ncols=2,
            label_box=SlideBox(0.35, 1.5, 2.5, 0.8),
            column_boxes=(
                SlideBox(0.35, 2.1, 1.1, 11.5),
                SlideBox(1.55, 2.1, 1.1, 11.5),
            ),
        ),
        SlideGroupSpec(
            label="3 Atoms",
            natoms=3,
            ncols=2,
            label_box=SlideBox(2.85, 1.5, 2.5, 0.8),
            column_boxes=(
                SlideBox(2.85, 2.1, 1.4, 11.0),
                SlideBox(4.15, 2.1, 1.4, 11.0),
            ),
        ),
        SlideGroupSpec(
            label="4 Atoms",
            natoms=4,
            ncols=2,
            label_box=SlideBox(5.85, 1.5, 2.5, 0.8),
            column_boxes=(
                SlideBox(5.85, 2.1, 1.6, 8.0),
                SlideBox(7.45, 2.1, 1.6, 8.0),
            ),
        ),
        SlideGroupSpec(
            label="5 Atoms",
            natoms=5,
            ncols=2,
            label_box=SlideBox(9.25, 1.5, 2.5, 0.8),
            column_boxes=(
                SlideBox(9.25, 2.1, 1.9, 8.0),
                SlideBox(11.15, 2.1, 1.9, 8.0),
            ),
        ),
        SlideGroupSpec(
            label="6 Atoms",
            natoms=6,
            ncols=1,
            label_box=SlideBox(13.25, 1.5, 2.5, 0.8),
            column_boxes=(SlideBox(13.25, 2.1, 2.5, 10.5),),
        ),
        SlideGroupSpec(
            label="7 Atoms",
            natoms=7,
            ncols=1,
            label_box=SlideBox(15.55, 1.5, 2.5, 0.8),
            column_boxes=(SlideBox(15.55, 2.1, 2.5, 10.5),),
        ),
        SlideGroupSpec(
            label="8 Atoms",
            natoms=8,
            ncols=1,
            label_box=SlideBox(18.05, 1.5, 2.5, 0.8),
            column_boxes=(SlideBox(18.05, 2.1, 2.5, 10.5),),
        ),
        SlideGroupSpec(
            label="9 Atoms",
            natoms=9,
            ncols=2,
            label_box=SlideBox(20.75, 1.5, 2.5, 0.8),
            column_boxes=(
                SlideBox(20.75, 2.1, 3.0, 4.0),
                SlideBox(23.25, 2.1, 3.0, 4.0),
            ),
        ),
        SlideGroupSpec(
            label="10 Atoms",
            natoms=10,
            ncols=1,
            label_box=SlideBox(20.75, 5.2, 2.5, 0.8),
            column_boxes=(SlideBox(20.75, 5.65, 3.0, 3.0),),
        ),
        SlideGroupSpec(
            label="11 Atoms",
            natoms=11,
            ncols=1,
            label_box=SlideBox(23.25, 5.2, 2.5, 0.8),
            column_boxes=(SlideBox(23.25, 5.65, 3.0, 3.0),),
        ),
        SlideGroupSpec(
            label="12 Atoms",
            natoms=12,
            ncols=1,
            label_box=SlideBox(15.55, 9.6, 2.5, 0.8),
            column_boxes=(SlideBox(15.55, 10.1, 3.0, 2.0),),
        ),
        SlideGroupSpec(
            label="13+ Atoms",
            natoms=13,
            ncols=2,
            label_box=SlideBox(20.75, 10.45, 2.5, 0.8),
            column_boxes=(
                SlideBox(20.75, 11.0, 3.0, 2.0),
                SlideBox(23.25, 11.0, 3.0, 2.0),
            ),
            natoms_greater=True,
        ),
    )


def balanced_molecule_slide_specs(
    entries: list[SlideMoleculeEntry],
) -> tuple[tuple[SlideGroupSpec, ...], float]:
    """Return adaptive slide group specs and the chosen molecule font size.

    The balanced profile keeps all atom-count groups in a single molecule band
    and adapts by choosing a readable font size, deriving column counts from
    row capacity, estimating per-group widths from display formulas, and
    distributing remaining horizontal space as group gaps.
    """
    for font_pt in BALANCED_MOLECULE_FONT_CANDIDATES:
        specs = _balanced_zoned_specs_for_font(entries, font_pt)
        if _balanced_zones_fit(specs):
            return specs, font_pt

    font_pt = BALANCED_MOLECULE_FONT_CANDIDATES[-1]
    return _balanced_zoned_specs_for_font(entries, font_pt), font_pt


def compact_molecule_slide_specs(
    entries: list[SlideMoleculeEntry],
) -> tuple[tuple[SlideGroupSpec, ...], float]:
    """Return compact slide specs for sparse molecule inventories."""
    for font_pt in COMPACT_MOLECULE_FONT_CANDIDATES:
        specs = _compact_specs_for_font(entries, font_pt)
        if _compact_specs_fit(specs):
            return specs, font_pt

    font_pt = COMPACT_MOLECULE_FONT_CANDIDATES[-1]
    return _compact_specs_for_font(entries, font_pt), font_pt


def _compact_specs_for_font(
    entries: list[SlideMoleculeEntry],
    font_pt: float,
) -> tuple[SlideGroupSpec, ...]:
    """Return specs for occupied atom-count bins in a single compact band."""
    groups = _entries_by_atom_group(entries)
    occupied_natoms = tuple(
        natoms
        for natoms in sorted(groups)
        if groups[natoms]
    )
    specs = _left_packed_compact_specs(
        groups,
        occupied_natoms=occupied_natoms,
        font_pt=font_pt,
    )
    return _spread_compact_specs(specs)


def _left_packed_compact_specs(
    groups: dict[int, list[SlideMoleculeEntry]],
    *,
    occupied_natoms: tuple[int, ...],
    font_pt: float,
) -> tuple[SlideGroupSpec, ...]:
    """Return left-packed compact specs for occupied atom-count bins."""
    column_height = COMPACT_COLUMN_BOTTOM - COMPACT_COLUMN_TOP
    capacity = max(
        1,
        estimated_column_capacity(
            SlideBox(0, 0, 1, column_height),
            font_pt=font_pt,
        ),
    )

    specs = []
    cursor = COMPACT_CONTENT_LEFT
    for natoms in occupied_natoms:
        label = f"{natoms} Atoms" if natoms < 13 else "13+ Atoms"
        group_entries = groups[natoms]
        ncols = max(1, ceil(len(group_entries) / capacity))
        column_chunks = _split_entries_for_columns(group_entries, ncols)
        column_widths = [
            balanced_column_width(chunk, font_pt)
            for chunk in column_chunks
        ]
        columns_width = sum(column_widths) + max(0, ncols - 1) * BALANCED_COLUMN_GAP
        group_width = max(columns_width, compact_label_width(label))
        offset = (group_width - columns_width) / 2
        column_left = cursor + offset
        column_boxes = []
        for width in column_widths:
            column_boxes.append(
                SlideBox(
                    column_left,
                    COMPACT_COLUMN_TOP,
                    width,
                    column_height,
                )
            )
            column_left += width + BALANCED_COLUMN_GAP
        specs.append(
            SlideGroupSpec(
                label=label,
                natoms=natoms,
                ncols=ncols,
                label_box=SlideBox(cursor, COMPACT_LABEL_TOP, group_width, 0.5),
                column_boxes=tuple(column_boxes),
                natoms_greater=natoms == 13,
            )
        )
        cursor += group_width + COMPACT_GROUP_GAP
    return tuple(specs)


def _spread_compact_specs(
    specs: tuple[SlideGroupSpec, ...],
) -> tuple[SlideGroupSpec, ...]:
    """Distribute compact groups across the slide without creating large gaps."""
    return _spread_zone_specs(
        specs,
        left=COMPACT_CONTENT_LEFT,
        right=COMPACT_CONTENT_RIGHT,
        base_gap=COMPACT_GROUP_GAP,
        max_gap=COMPACT_MAX_GROUP_GAP,
    )


def _compact_specs_fit(specs: tuple[SlideGroupSpec, ...]) -> bool:
    """Return whether compact specs fit inside the compact slide band."""
    if not specs:
        return True
    return (
        min(spec.label_box.left for spec in specs) >= COMPACT_CONTENT_LEFT
        and max(spec.label_box.right for spec in specs) <= COMPACT_CONTENT_RIGHT
    )


def _balanced_zoned_specs_for_font(
    entries: list[SlideMoleculeEntry],
    font_pt: float,
) -> tuple[SlideGroupSpec, ...]:
    """Return balanced specs using a top row and bottom-right large-molecule row."""
    groups = _entries_by_atom_group(entries)
    top_specs = _balanced_zone_specs(
        groups,
        natoms_values=BALANCED_TOP_NATOMS,
        font_pt=font_pt,
        label_top=BALANCED_TOP_LABEL_TOP,
        column_top=BALANCED_TOP_COLUMN_TOP,
        column_bottom=BALANCED_TOP_COLUMN_BOTTOM,
        left=BALANCED_CONTENT_LEFT,
        right=BALANCED_CONTENT_RIGHT,
        min_columns=BALANCED_TOP_MIN_COLUMNS,
    )
    bottom_specs = _balanced_zone_specs(
        groups,
        natoms_values=BALANCED_BOTTOM_NATOMS,
        font_pt=font_pt,
        label_top=BALANCED_BOTTOM_LABEL_TOP,
        column_top=BALANCED_BOTTOM_COLUMN_TOP,
        column_bottom=BALANCED_BOTTOM_COLUMN_BOTTOM,
        left=BALANCED_BOTTOM_LEFT,
        right=BALANCED_CONTENT_RIGHT,
        min_columns={},
    )
    bottom_specs = _align_bottom_specs_to_top_columns(top_specs, bottom_specs)
    return top_specs + bottom_specs


def _balanced_zone_specs(
    groups: dict[int, list[SlideMoleculeEntry]],
    *,
    natoms_values: tuple[int, ...],
    font_pt: float,
    label_top: float,
    column_top: float,
    column_bottom: float,
    left: float,
    right: float,
    min_columns: dict[int, int],
) -> tuple[SlideGroupSpec, ...]:
    """Return spread group specs for one horizontal layout zone."""
    specs = _left_packed_balanced_zone_specs(
        groups,
        natoms_values=natoms_values,
        font_pt=font_pt,
        label_top=label_top,
        column_top=column_top,
        column_bottom=column_bottom,
        left=left,
        min_columns=min_columns,
    )
    return _spread_zone_specs(specs, left=left, right=right)


def _left_packed_balanced_zone_specs(
    groups: dict[int, list[SlideMoleculeEntry]],
    *,
    natoms_values: tuple[int, ...],
    font_pt: float,
    label_top: float,
    column_top: float,
    column_bottom: float,
    left: float,
    min_columns: dict[int, int],
) -> tuple[SlideGroupSpec, ...]:
    """Return left-packed specs for one balanced layout zone."""
    column_height = column_bottom - column_top
    capacity = max(
        1,
        estimated_column_capacity(
            SlideBox(0, 0, 1, column_height),
            font_pt=font_pt,
        ),
    )

    specs = []
    cursor = left
    for natoms in natoms_values:
        label = f"{natoms} Atoms" if natoms < 13 else "13+ Atoms"
        group_entries = groups[natoms]
        ncols = max(
            min_columns.get(natoms, 1),
            ceil(len(group_entries) / capacity) if group_entries else 1,
        )
        column_chunks = _split_entries_for_columns(group_entries, ncols)
        column_widths = [
            balanced_column_width(chunk, font_pt)
            for chunk in column_chunks
        ]
        columns_width = sum(column_widths) + max(0, ncols - 1) * BALANCED_COLUMN_GAP
        group_width = max(columns_width, balanced_label_width(label))
        offset = (group_width - columns_width) / 2
        column_left = cursor + offset
        column_boxes = []
        for width in column_widths:
            column_boxes.append(
                SlideBox(
                    column_left,
                    column_top,
                    width,
                    column_height,
                )
            )
            column_left += width + BALANCED_COLUMN_GAP
        column_boxes = tuple(
            column_boxes
        )
        specs.append(
            SlideGroupSpec(
                label=label,
                natoms=natoms,
                ncols=ncols,
                label_box=SlideBox(cursor, label_top, group_width, 0.45),
                column_boxes=column_boxes,
                natoms_greater=natoms == 13,
            )
        )
        cursor += group_width + BALANCED_GROUP_GAP
    return tuple(specs)


def _align_bottom_specs_to_top_columns(
    top_specs: tuple[SlideGroupSpec, ...],
    bottom_specs: tuple[SlideGroupSpec, ...],
) -> tuple[SlideGroupSpec, ...]:
    """Align the large-molecule row with selected upper-row column starts."""
    top_by_natoms = {spec.natoms: spec for spec in top_specs}
    anchors = {
        11: (top_by_natoms[6].column_boxes[0].left,),
        12: (top_by_natoms[7].column_boxes[0].left,),
        13: (
            top_by_natoms[9].column_boxes[0].left,
            top_by_natoms[10].column_boxes[0].left,
        ),
    }

    aligned = []
    for spec in bottom_specs:
        if spec.natoms not in anchors:
            aligned.append(spec)
            continue

        anchor_values = anchors[spec.natoms]
        column_boxes = []
        for index, box in enumerate(spec.column_boxes):
            left = anchor_values[min(index, len(anchor_values) - 1)]
            column_boxes.append(
                SlideBox(
                    left,
                    box.top,
                    box.width,
                    box.height,
                )
            )

        label_left = column_boxes[0].left
        label_right = max(box.right for box in column_boxes)
        label_width = max(
            label_right - label_left,
            balanced_label_width(spec.label),
        )
        aligned.append(
            SlideGroupSpec(
                label=spec.label,
                natoms=spec.natoms,
                ncols=spec.ncols,
                label_box=SlideBox(
                    label_left,
                    spec.label_box.top,
                    label_width,
                    spec.label_box.height,
                ),
                column_boxes=tuple(column_boxes),
                natoms_greater=spec.natoms_greater,
            )
        )
    return tuple(aligned)


def balanced_total_box(specs: tuple[SlideGroupSpec, ...]) -> SlideBox:
    """Return a count box centered under the first six molecule columns."""
    box_width = 5.9
    anchor_natoms = {2, 3, 4}
    columns = [
        column
        for spec in specs
        if spec.natoms in anchor_natoms
        for column in spec.column_boxes
    ]
    if not columns:
        return SlideBox(0.8, 10.75, 5.6, 1.35)

    left = min(column.left for column in columns)
    right = max(column.right for column in columns)
    center = (left + right) / 2
    return SlideBox(center - box_width / 2, 10.75, box_width, 1.28)


def compact_total_box(specs: tuple[SlideGroupSpec, ...]) -> SlideBox:
    """Return a centered count box for a compact slide layout."""
    box_width = 5.9
    if not specs:
        return SlideBox((LEGACY_SLIDE_WIDTH_IN - box_width) / 2, 12.05, box_width, 1.28)

    left = min(column.left for spec in specs for column in spec.column_boxes)
    right = max(column.right for spec in specs for column in spec.column_boxes)
    center = (left + right) / 2
    return SlideBox(center - box_width / 2, 12.05, box_width, 1.28)


def rendered_group_label_box(group: SlideGroupPlan, profile: str) -> SlideBox:
    """Return the textbox used for rendering a group label."""
    if profile not in {"balanced", "compact"}:
        return group.spec.label_box

    left = min(column.box.left for column in group.columns)
    right = max(column.box.right for column in group.columns)
    label_width = (
        compact_label_width(group.spec.label)
        if profile == "compact"
        else balanced_label_width(group.spec.label)
    )
    if len(group.columns) == 1:
        width = max(right - left, label_width)
        return SlideBox(
            left,
            group.spec.label_box.top,
            width,
            group.spec.label_box.height,
        )

    return SlideBox(
        left,
        group.spec.label_box.top,
        right - left,
        group.spec.label_box.height,
    )


def group_label_alignment(group: SlideGroupPlan, profile: str, PP_ALIGN):
    """Return label alignment for a group header."""
    if profile in {"balanced", "compact"} and len(group.columns) == 1:
        return PP_ALIGN.LEFT
    if profile in {"balanced", "compact"} and len(group.columns) > 1:
        return PP_ALIGN.CENTER
    return None


def _split_entries_for_columns(
    entries: list[SlideMoleculeEntry],
    ncols: int,
) -> list[list[SlideMoleculeEntry]]:
    """Split entries using the same rule as the rendered group planner."""
    if ncols <= 0:
        return []
    chunk_size = ceil(len(entries) / ncols) if entries else 0
    return [
        entries[index * chunk_size : (index + 1) * chunk_size]
        for index in range(ncols)
    ]


def _spread_zone_specs(
    specs: tuple[SlideGroupSpec, ...],
    *,
    left: float,
    right: float,
    base_gap: float = BALANCED_GROUP_GAP,
    max_gap: float = BALANCED_MAX_GROUP_GAP,
) -> tuple[SlideGroupSpec, ...]:
    """Distribute extra horizontal space within a layout zone."""
    if not specs:
        return specs

    used = max(spec.label_box.right for spec in specs) - left
    available = right - left
    if len(specs) == 1:
        dx = max(0, available - used) / 2
        return (_move_group_spec(specs[0], dx),)
    if used >= available:
        return specs

    extra = available - used
    extra_per_gap = min(
        extra / (len(specs) - 1),
        max_gap - base_gap,
    )
    leftover = extra - extra_per_gap * (len(specs) - 1)
    cursor = left + leftover / 2
    spread = []
    for spec in specs:
        dx = cursor - spec.label_box.left
        spread.append(_move_group_spec(spec, dx))
        cursor += spec.label_box.width + base_gap + extra_per_gap
    return tuple(spread)


def _balanced_zones_fit(specs: tuple[SlideGroupSpec, ...]) -> bool:
    """Return whether all balanced zones are within their intended bounds."""
    top_specs = [spec for spec in specs if spec.natoms in BALANCED_TOP_NATOMS]
    bottom_specs = [spec for spec in specs if spec.natoms in BALANCED_BOTTOM_NATOMS]
    if not top_specs or not bottom_specs:
        return False
    top_left = min(spec.label_box.left for spec in top_specs)
    top_right = max(spec.label_box.right for spec in top_specs)
    bottom_left = min(spec.label_box.left for spec in bottom_specs)
    bottom_right = max(spec.label_box.right for spec in bottom_specs)
    return (
        top_left >= BALANCED_CONTENT_LEFT
        and top_right <= BALANCED_CONTENT_RIGHT
        and bottom_left >= BALANCED_CONTENT_LEFT
        and bottom_right <= BALANCED_CONTENT_RIGHT
    )


def _balanced_specs_for_font(
    entries: list[SlideMoleculeEntry],
    font_pt: float,
) -> tuple[SlideGroupSpec, ...]:
    """Return left-packed balanced specs for one candidate font size."""
    groups = _entries_by_atom_group(entries)
    column_height = BALANCED_COLUMN_BOTTOM - BALANCED_COLUMN_TOP
    capacity = max(
        1,
        estimated_column_capacity(
            SlideBox(0, 0, 1, column_height),
            font_pt=font_pt,
        ),
    )

    specs = []
    left = BALANCED_CONTENT_LEFT
    for natoms in list(range(2, 13)) + [13]:
        label = f"{natoms} Atoms" if natoms < 13 else "13+ Atoms"
        group_entries = groups[natoms]
        ncols = max(1, ceil(len(group_entries) / capacity))
        column_width = balanced_column_width(group_entries, font_pt)
        column_boxes = tuple(
            SlideBox(
                left + index * (column_width + BALANCED_COLUMN_GAP),
                BALANCED_COLUMN_TOP,
                column_width,
                column_height,
            )
            for index in range(ncols)
        )
        columns_width = (
            ncols * column_width
            + max(0, ncols - 1) * BALANCED_COLUMN_GAP
        )
        group_width = max(columns_width, balanced_label_width(label))
        offset = (group_width - columns_width) / 2
        column_boxes = tuple(
            SlideBox(
                box.left + offset,
                box.top,
                box.width,
                box.height,
            )
            for box in column_boxes
        )
        specs.append(
            SlideGroupSpec(
                label=label,
                natoms=natoms,
                ncols=ncols,
                label_box=SlideBox(left, BALANCED_LABEL_TOP, group_width, 0.45),
                column_boxes=column_boxes,
                natoms_greater=natoms == 13,
            )
        )
        left += group_width + BALANCED_GROUP_GAP
    return tuple(specs)


def _spread_balanced_specs(
    specs: tuple[SlideGroupSpec, ...],
) -> tuple[SlideGroupSpec, ...]:
    """Distribute extra horizontal space between balanced groups."""
    available = _balanced_available_width()
    used = _balanced_specs_width(specs)
    if len(specs) <= 1 or used >= available:
        return specs

    base_left = min(spec.label_box.left for spec in specs)
    extra = available - used
    extra_per_gap = min(extra / (len(specs) - 1), BALANCED_MAX_GROUP_GAP - BALANCED_GROUP_GAP)
    leftover = extra - extra_per_gap * (len(specs) - 1)
    left = base_left + leftover / 2

    spread = []
    for spec in specs:
        dx = left - spec.label_box.left
        spread.append(_move_group_spec(spec, dx))
        left += spec.label_box.width + BALANCED_GROUP_GAP + extra_per_gap
    return tuple(spread)


def _balanced_specs_width(specs: tuple[SlideGroupSpec, ...]) -> float:
    """Return occupied balanced-spec width including current gaps."""
    if not specs:
        return 0.0
    left = min(spec.label_box.left for spec in specs)
    right = max(spec.label_box.right for spec in specs)
    return right - left


def _balanced_available_width() -> float:
    """Return horizontal width available to the balanced molecule band."""
    return BALANCED_CONTENT_RIGHT - BALANCED_CONTENT_LEFT


def _move_group_spec(spec: SlideGroupSpec, dx: float) -> SlideGroupSpec:
    """Move a group spec horizontally."""
    return SlideGroupSpec(
        label=spec.label,
        natoms=spec.natoms,
        ncols=spec.ncols,
        label_box=SlideBox(
            spec.label_box.left + dx,
            spec.label_box.top,
            spec.label_box.width,
            spec.label_box.height,
        ),
        column_boxes=tuple(
            SlideBox(
                box.left + dx,
                box.top,
                box.width,
                box.height,
            )
            for box in spec.column_boxes
        ),
        natoms_greater=spec.natoms_greater,
    )


def _entries_by_atom_group(
    entries: list[SlideMoleculeEntry],
) -> dict[int, list[SlideMoleculeEntry]]:
    """Group slide entries into 2-12 atom bins and a 13+ bin keyed as 13."""
    groups = {natoms: [] for natoms in range(2, 14)}
    for entry in entries:
        key = entry.natoms if entry.natoms < 13 else 13
        groups.setdefault(key, []).append(entry)
    return groups


def balanced_column_width(
    entries: list[SlideMoleculeEntry],
    font_pt: float,
) -> float:
    """Estimate a readable PowerPoint column width for formula entries."""
    visible_lengths = [
        len("".join(run.text for run in powerpoint_formula_runs(entry.display_formula)))
        for entry in entries
    ]
    max_length = max(visible_lengths, default=5)
    return max(
        0.95,
        min(3.0, 0.22 + max_length * font_pt * BALANCED_FORMULA_WIDTH_FACTOR),
    )


def balanced_label_width(label: str) -> float:
    """Return a minimum group width that leaves room for the atom-count label."""
    return max(1.45, len(label) * 22 * 0.0075)


def compact_label_width(label: str) -> float:
    """Return a minimum compact-profile group width for one-line labels."""
    return max(1.85, len(label) * 24 * 0.010)


def selected_molecule_slide_entries(
    view: CensusView,
    *,
    detection_type: str = DEFAULT_DETECTION_TYPE,
    include_tentative: bool = False,
    include_disputed: bool = False,
    include_isotopologues: bool = False,
) -> list[SlideMoleculeEntry]:
    """Return molecules selected for a molecule-summary slide.

    The default is the main secure, non-isotopologue ISM/CSM inventory.
    Ordering follows first accepted detection order within the selected view,
    which is the order used for the manuscript molecule tables.
    """
    molecules_by_label = {
        molecule.label: molecule
        for molecule in view.context_molecules(
            detection_type,
            include_tentative=include_tentative,
            include_disputed=include_disputed,
            include_isotopologues=include_isotopologues,
        )
    }
    detections = view.context_detections(
        detection_type,
        include_tentative=include_tentative,
        include_disputed=include_disputed,
        include_isotopologues=include_isotopologues,
    )
    sort_keys = first_detection_sort_keys(detections)
    molecules = sorted(
        molecules_by_label.values(),
        key=lambda molecule: sort_keys.get(
            molecule.label,
            (date.max, molecule.natoms, molecule.label, ""),
        ),
    )
    return [
        SlideMoleculeEntry(
            molecule=molecule,
            display_formula=molecule.table_formula or molecule.formula,
        )
        for molecule in molecules
    ]


def first_detection_sort_keys(detections: list[Detection]) -> dict[str, tuple]:
    """Return first-detection sort keys by molecule label."""
    keys = {}
    for detection in detections:
        label = detection.molecule.label
        key = (
            detection.sortdate,
            detection.molecule.natoms,
            detection.molecule.label,
            detection.id,
        )
        if label not in keys or key < keys[label]:
            keys[label] = key
    return keys


def build_molecule_slide_layout(
    view: CensusView,
    *,
    detection_type: str = DEFAULT_DETECTION_TYPE,
    include_tentative: bool = False,
    include_disputed: bool = False,
    include_isotopologues: bool = False,
    profile: str = "legacy",
    title: str = LEGACY_TITLE,
    last_updated: str | None = None,
) -> MoleculeSlideLayout:
    """Build a molecule-slide layout plan without writing PowerPoint output."""
    entries = selected_molecule_slide_entries(
        view,
        detection_type=detection_type,
        include_tentative=include_tentative,
        include_disputed=include_disputed,
        include_isotopologues=include_isotopologues,
    )

    if profile == "legacy":
        molecule_font_pt = LEGACY_MOLECULE_FONT_PT
        specs = legacy_molecule_slide_specs()
        title_box = SlideBox(0.35, 0.22, 18.33, 1.18)
        credit_box = SlideBox(22.15, 0.22, 300 / 72, 90 / 72)
        total_box = SlideBox(6.35, 10.6, 6.0, 1.4)
        footer_box = SlideBox(6.35, 12.2, 6.0, 1.0)
        total_font_pt = 45
        footer_font_pt = 32
        total_line_width_pt = 6
    elif profile == "balanced":
        specs, molecule_font_pt = balanced_molecule_slide_specs(entries)
        title_box = SlideBox(0.35, 0.22, 18.33, 1.18)
        credit_box = SlideBox(22.15, 0.22, 300 / 72, 90 / 72)
        total_box = balanced_total_box(specs)
        footer_box = SlideBox(
            total_box.left,
            total_box.bottom + 0.08,
            total_box.width,
            0.5,
        )
        group_label_font_pt = 22
        total_font_pt = 43
        footer_font_pt = 25
        total_line_width_pt = 5
    elif profile == "compact":
        specs, molecule_font_pt = compact_molecule_slide_specs(entries)
        title_box = SlideBox(0.35, 0.22, 18.33, 1.18)
        credit_box = SlideBox(22.15, 0.22, 300 / 72, 90 / 72)
        total_box = compact_total_box(specs)
        footer_box = SlideBox(
            total_box.left,
            total_box.bottom + 0.08,
            total_box.width,
            0.5,
        )
        group_label_font_pt = 24
        total_font_pt = 43
        footer_font_pt = 25
        total_line_width_pt = 5
    else:
        raise ValueError(
            "Molecule slide profile must be 'legacy', 'balanced', or 'compact'."
        )

    groups = tuple(
        _group_plan(spec, entries, font_pt=molecule_font_pt)
        for spec in specs
    )
    warnings = tuple(_layout_warnings(groups, font_pt=molecule_font_pt))
    return MoleculeSlideLayout(
        groups=groups,
        title=title,
        total=len(entries),
        profile=profile,
        title_box=title_box,
        credit_box=credit_box,
        total_box=total_box,
        footer_box=footer_box,
        molecule_font_pt=molecule_font_pt,
        group_label_font_pt=group_label_font_pt if profile == "balanced" else 26,
        total_font_pt=total_font_pt,
        footer_font_pt=footer_font_pt,
        total_line_width_pt=total_line_width_pt,
        package_version=package_version(),
        last_updated=last_updated or selected_records_last_modified(entries),
        detection_type=detection_type,
        include_tentative=include_tentative,
        include_disputed=include_disputed,
        include_isotopologues=include_isotopologues,
        warnings=warnings,
    )


def _group_plan(
    spec: SlideGroupSpec,
    entries: list[SlideMoleculeEntry],
    *,
    font_pt: float = LEGACY_MOLECULE_FONT_PT,
) -> SlideGroupPlan:
    """Return a rendered group plan for a static group specification."""
    if spec.natoms_greater:
        group_entries = [entry for entry in entries if entry.natoms >= spec.natoms]
    else:
        group_entries = [entry for entry in entries if entry.natoms == spec.natoms]

    chunk_size = ceil(len(group_entries) / spec.ncols) if group_entries else 0
    columns = []
    for index in range(spec.ncols):
        start = index * chunk_size
        stop = start + chunk_size
        molecules = tuple(group_entries[start:stop])
        box = spec.column_boxes[index]
        columns.append(
            SlideColumnPlan(
                box=box,
                molecules=molecules,
                estimated_capacity=estimated_column_capacity(box, font_pt=font_pt),
            )
        )
    return SlideGroupPlan(spec=spec, columns=tuple(columns))


def estimated_column_capacity(
    box: SlideBox,
    *,
    font_pt: float = LEGACY_MOLECULE_FONT_PT,
) -> int:
    """Return a conservative estimated row capacity for a column box."""
    row_height_in = font_pt * LEGACY_ROW_SPACING / 72
    return int(box.height / row_height_in)


def _layout_warnings(
    groups: tuple[SlideGroupPlan, ...],
    *,
    font_pt: float = LEGACY_MOLECULE_FONT_PT,
) -> list[str]:
    """Return human-readable layout warnings."""
    warnings = []
    for group in groups:
        if group.spec.label_box.right > LEGACY_SLIDE_WIDTH_IN:
            warnings.append(f"{group.spec.label} label exceeds slide width.")
        if group.spec.label_box.bottom > LEGACY_SLIDE_HEIGHT_IN:
            warnings.append(f"{group.spec.label} label exceeds slide height.")
        for index, column in enumerate(group.columns, start=1):
            if column.overflow:
                warnings.append(
                    f"{group.spec.label} column {index} exceeds estimated "
                    f"capacity by {column.overflow} row(s)."
                )
            if column.box.right > LEGACY_SLIDE_WIDTH_IN:
                warnings.append(f"{group.spec.label} column {index} exceeds slide width.")
            if column.box.bottom > LEGACY_SLIDE_HEIGHT_IN:
                warnings.append(f"{group.spec.label} column {index} exceeds slide height.")
    return warnings


def selected_records_last_modified(entries: list[SlideMoleculeEntry]) -> str | None:
    """Return the latest molecule-history modification date in selected entries."""
    values = [
        entry.molecule.history.last_modified
        for entry in entries
        if entry.molecule.history is not None
        and entry.molecule.history.last_modified is not None
    ]
    return max(values) if values else None


def package_version() -> str:
    """Return installed package version, or a traceable development label."""
    try:
        return metadata.version("astromol")
    except metadata.PackageNotFoundError:
        git_label = repository_git_label()
        if git_label is None:
            return "development"
        return f"development (git {git_label})"


def repository_git_label() -> str | None:
    """Return the current repository short hash, with dirty marker if needed."""
    repo_root = Path(__file__).resolve().parents[1]
    try:
        head = subprocess.run(
            ["git", "-C", str(repo_root), "rev-parse", "--short", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
        status = subprocess.run(
            ["git", "-C", str(repo_root), "status", "--porcelain"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
    except (FileNotFoundError, subprocess.SubprocessError):
        return None
    if not head:
        return None
    return f"{head}+dirty" if status else head


def slide_version_label(version: str) -> str:
    """Return display text for the slide package-version line."""
    if version.startswith("development"):
        return version
    return f"v{version}"


def molecule_slide_report(layout: MoleculeSlideLayout) -> str:
    """Return a concise Markdown report for a molecule slide layout."""
    lines = [
        "# Molecule Slide Layout Report",
        "",
        f"- profile: `{layout.profile}`",
        f"- detection type: `{layout.detection_type}`",
        f"- total molecules: {layout.total}",
        f"- molecule font: {layout.molecule_font_pt:g} pt",
        f"- include tentative: {layout.include_tentative}",
        f"- include disputed: {layout.include_disputed}",
        f"- include isotopologues: {layout.include_isotopologues}",
        f"- last updated: {layout.last_updated or 'not available'}",
        "",
        "## Groups",
        "",
    ]
    for group in layout.groups:
        column_counts = ", ".join(str(len(column.molecules)) for column in group.columns)
        lines.append(
            f"- {group.spec.label}: {len(group.molecules)} molecule(s); "
            f"columns = {column_counts}"
        )
    lines.extend(["", "## Warnings", ""])
    if layout.warnings:
        lines.extend(f"- {warning}" for warning in layout.warnings)
    else:
        lines.append("- none")
    return "\n".join(lines) + "\n"


def write_molecule_slide_report(
    layout: MoleculeSlideLayout,
    output_path: str | Path,
) -> Path:
    """Write a molecule-slide layout report and return its path."""
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(molecule_slide_report(layout))
    return path


def write_molecule_slide(
    view: CensusView,
    output_path: str | Path = "astro_molecules.pptx",
    *,
    detection_type: str = DEFAULT_DETECTION_TYPE,
    include_tentative: bool = False,
    include_disputed: bool = False,
    include_isotopologues: bool = False,
    profile: str = "legacy",
    title: str = LEGACY_TITLE,
    last_updated: str | None = None,
) -> Path:
    """Write a molecule-summary PowerPoint slide and return the output path."""
    layout = build_molecule_slide_layout(
        view,
        detection_type=detection_type,
        include_tentative=include_tentative,
        include_disputed=include_disputed,
        include_isotopologues=include_isotopologues,
        profile=profile,
        title=title,
        last_updated=last_updated,
    )
    presentation = render_molecule_slide(layout)
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    presentation.save(path)
    return path


def build_ppd_detection_slide_layout(
    view: CensusView,
    *,
    include_tentative: bool = False,
    include_disputed: bool = False,
    include_isotopologues: bool = True,
    profile: str = "compact",
    title: str = PPD_TITLE,
    last_updated: str | None = None,
) -> MoleculeSlideLayout:
    """Build the PPD molecule/isotopologue slide layout plan."""
    return build_molecule_slide_layout(
        view,
        detection_type=PPD_DETECTION_TYPE,
        include_tentative=include_tentative,
        include_disputed=include_disputed,
        include_isotopologues=include_isotopologues,
        profile=profile,
        title=title,
        last_updated=last_updated,
    )


def write_ppd_detection_slide(
    view: CensusView,
    output_path: str | Path = "ppd_molecules.pptx",
    *,
    include_tentative: bool = False,
    include_disputed: bool = False,
    include_isotopologues: bool = True,
    profile: str = "compact",
    title: str = PPD_TITLE,
    last_updated: str | None = None,
) -> Path:
    """Write a PPD molecule/isotopologue PowerPoint slide."""
    layout = build_ppd_detection_slide_layout(
        view,
        include_tentative=include_tentative,
        include_disputed=include_disputed,
        include_isotopologues=include_isotopologues,
        profile=profile,
        title=title,
        last_updated=last_updated,
    )
    presentation = render_molecule_slide(layout)
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    presentation.save(path)
    return path


def render_molecule_slide(layout: MoleculeSlideLayout):
    """Render a molecule-slide layout to a ``python-pptx`` Presentation."""
    try:
        from pptx import Presentation
        from pptx.dml.color import RGBColor
        from pptx.enum.shapes import MSO_SHAPE
        from pptx.enum.text import MSO_ANCHOR, PP_ALIGN
        from pptx.util import Inches, Pt
    except ImportError as exc:
        raise ImportError(
            "PowerPoint slide generation requires the 'python-pptx' package."
        ) from exc

    presentation = Presentation()
    presentation.slide_width = Pt(layout.slide_width_pt)
    presentation.slide_height = Pt(layout.slide_height_pt)

    slide = presentation.slides.add_slide(presentation.slide_layouts[6])
    _add_text_box(
        slide,
        Inches,
        Pt,
        layout.title_box,
        layout.title,
        font_size=64,
        bold=False,
    )
    _add_credit_box(slide, Inches, Pt, PP_ALIGN, layout)

    for group in layout.groups:
        label_box = rendered_group_label_box(group, layout.profile)
        _add_text_box(
            slide,
            Inches,
            Pt,
            label_box,
            group.spec.label,
            font_size=layout.group_label_font_pt,
            bold=True,
            align=group_label_alignment(group, layout.profile, PP_ALIGN),
        )
        for column in group.columns:
            _add_molecule_column(slide, Inches, Pt, column, layout.molecule_font_pt)

    _add_total_box(
        slide,
        Inches,
        Pt,
        RGBColor,
        MSO_SHAPE,
        MSO_ANCHOR,
        PP_ALIGN,
        layout,
    )
    if layout.last_updated:
        _add_text_box(
            slide,
            Inches,
            Pt,
            layout.footer_box,
            f"Last Updated: {format_slide_date(layout.last_updated)}",
            font_size=layout.footer_font_pt,
            bold=True,
            align=PP_ALIGN.CENTER,
        )
    return presentation


def _add_text_box(
    slide,
    Inches,
    Pt,
    box: SlideBox,
    text: str,
    *,
    font_size: int,
    bold: bool = False,
    align=None,
):
    """Add a single-run textbox to a slide."""
    shape = slide.shapes.add_textbox(
        Inches(box.left),
        Inches(box.top),
        Inches(box.width),
        Inches(box.height),
    )
    paragraph = shape.text_frame.paragraphs[0]
    paragraph.text = text
    paragraph.font.name = "Arial"
    paragraph.font.size = Pt(font_size)
    paragraph.font.bold = bold
    if align is not None:
        paragraph.alignment = align
    return shape


def _add_credit_box(slide, Inches, Pt, PP_ALIGN, layout: MoleculeSlideLayout) -> None:
    """Add the top-right ASTROMOL credit box."""
    shape = slide.shapes.add_textbox(
        Inches(layout.credit_box.left),
        Inches(layout.credit_box.top),
        Inches(layout.credit_box.width),
        Inches(layout.credit_box.height),
    )
    text_frame = shape.text_frame
    text_frame.clear()

    lines = [
        (f"Created with ASTROMOL {slide_version_label(layout.package_version)}", 16),
        (LEGACY_URL, 16),
        (LEGACY_CITATION, 18),
    ]
    for index, (text, size) in enumerate(lines):
        paragraph = text_frame.paragraphs[0] if index == 0 else text_frame.add_paragraph()
        paragraph.text = text
        paragraph.alignment = PP_ALIGN.RIGHT
        paragraph.font.name = "Arial"
        paragraph.font.size = Pt(size)
        if text == LEGACY_CITATION:
            paragraph.font.italic = True


def _add_molecule_column(
    slide,
    Inches,
    Pt,
    column: SlideColumnPlan,
    font_pt: float,
) -> None:
    """Add one column of formatted molecule formulas."""
    shape = slide.shapes.add_textbox(
        Inches(column.box.left),
        Inches(column.box.top),
        Inches(column.box.width),
        Inches(column.box.height),
    )
    text_frame = shape.text_frame
    text_frame.clear()
    text_frame.word_wrap = False
    text_frame.margin_left = Pt(0)
    text_frame.margin_right = Pt(0)
    for index, molecule in enumerate(column.molecules):
        paragraph = text_frame.paragraphs[0] if index == 0 else text_frame.add_paragraph()
        paragraph.font.name = "Arial"
        paragraph.font.size = Pt(font_pt)
        for run_spec in powerpoint_formula_runs(molecule.display_formula):
            run = paragraph.add_run()
            run.text = run_spec.text
            run.font.name = "Arial"
            run.font.size = Pt(font_pt)
            run.font.italic = run_spec.italic
            if run_spec.baseline == "subscript":
                run.font._element.set("baseline", "-25000")
            elif run_spec.baseline == "superscript":
                run.font._element.set("baseline", "30000")


def _add_total_box(
    slide,
    Inches,
    Pt,
    RGBColor,
    MSO_SHAPE,
    MSO_ANCHOR,
    PP_ALIGN,
    layout: MoleculeSlideLayout,
) -> None:
    """Add the historical molecule-count box."""
    shape = slide.shapes.add_shape(
        MSO_SHAPE.ROUNDED_RECTANGLE,
        Inches(layout.total_box.left),
        Inches(layout.total_box.top),
        Inches(layout.total_box.width),
        Inches(layout.total_box.height),
    )
    shape.fill.solid()
    shape.fill.fore_color.rgb = RGBColor(194, 192, 191)
    shape.line.color.rgb = RGBColor(0, 0, 0)
    shape.line.width = Pt(layout.total_line_width_pt)

    text_frame = shape.text_frame
    text_frame.vertical_anchor = MSO_ANCHOR.MIDDLE
    paragraph = text_frame.paragraphs[0]
    paragraph.text = f"{layout.total} Molecules"
    paragraph.alignment = PP_ALIGN.CENTER
    paragraph.font.name = "Arial"
    paragraph.font.size = Pt(layout.total_font_pt)
    paragraph.font.bold = True
    paragraph.font.color.rgb = RGBColor(163, 31, 52)


def format_slide_date(value: str | None) -> str:
    """Format ISO dates for slide display."""
    if value is None:
        return ""
    try:
        parsed = datetime.strptime(value, "%Y-%m-%d").date()
    except ValueError:
        return value
    return f"{parsed.day} {parsed.strftime('%B')} {parsed.year}"


def powerpoint_formula_runs(formula: str) -> list[FormulaRun]:
    """Return formatted runs for a molecule formula in PowerPoint text.

    This is not a full LaTeX or mhchem parser. It handles the formula forms used
    by astromol table formulas: atom counts, leading isomer/conformer labels,
    terminal charges, and isotope forms such as ``^{13}`` or ``[13C]``.
    """
    formula = strip_formula_wrappers(formula)
    tokens = [token for token in FORMULA_TOKEN.findall(formula) if token not in {"{", "}"}]
    runs: list[FormulaRun] = []

    for index, token in enumerate(tokens):
        if token in {r"\ce{", r"\mathrm{", r"\text{"}:
            continue
        if token.startswith("^{") and token.endswith("}"):
            runs.append(FormulaRun(token[2:-1], baseline="superscript"))
            continue
        if token.startswith("[") and token.endswith("]"):
            match = re.fullmatch(r"\[(\d+)([A-Z][a-z]?)\]", token)
            if match:
                runs.append(FormulaRun(match.group(1), baseline="superscript"))
                runs.append(FormulaRun(match.group(2)))
                continue
        if token.isdigit() and is_atom_count_token(tokens, index):
            runs.append(FormulaRun(token, baseline="subscript"))
            continue
        if token in {"+", "-"} and is_terminal_charge(tokens, index):
            runs.append(FormulaRun(token, baseline="superscript"))
            continue

        italic = index == 0 and token.islower()
        runs.append(FormulaRun(token, italic=italic))

    return runs


def strip_formula_wrappers(formula: str) -> str:
    """Remove common formula wrappers not meaningful in PowerPoint text."""
    text = formula.strip()
    for prefix in (r"\ce{", r"\mathrm{", r"\text{"):
        if text.startswith(prefix) and text.endswith("}"):
            return text[len(prefix) : -1]
    return text


def is_atom_count_token(tokens: list[str], index: int) -> bool:
    """Return whether a numeric token should be shown as an atom subscript."""
    if index == 0:
        return False
    previous = tokens[index - 1]
    if previous in {"+", "-", "^{", "{", "}"}:
        return False
    if previous.startswith("^{"):
        return False
    if index + 1 < len(tokens) and tokens[index + 1] == "-":
        return False
    return bool(re.fullmatch(r"[A-Za-z]+", previous))


def is_terminal_charge(tokens: list[str], index: int) -> bool:
    """Return whether a plus/minus token is a terminal molecular charge."""
    return index == len(tokens) - 1
