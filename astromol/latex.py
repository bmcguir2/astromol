"""LaTeX output helpers for census manuscripts."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from datetime import date
from math import ceil
from pathlib import Path

import numpy as np

from .census import CensusView
from .models import Detection, Molecule


RADIO_WAVELENGTHS = {"cm", "mm", "sub-mm"}
MOLREF_COMMAND = r"\providecommand{\molref}[2]{\hyperref[#1]{\ce{#2}}}"
BALANCED_ISM_MAX_COLUMNS = 7
BALANCED_ISM_MAX_ROWS = 23


ISM_TABLE_TWO_SEVEN_SPEC = (
    r"\begin{tabular*}{\textwidth}{l l @{\extracolsep{\fill}} l l  "
    r"@{\extracolsep{\fill}} l l  @{\extracolsep{\fill}} l l "
    r"@{\extracolsep{\fill}} l @{\extracolsep{\fill}} l}"
)

ISM_TABLE_EIGHT_MORE_SPEC = (
    r"\begin{tabular*}{\textwidth}{l @{\extracolsep{\fill}} l "
    r"@{\extracolsep{\fill}} l @{\extracolsep{\fill}} l "
    r"@{\extracolsep{\fill}} l @{\extracolsep{\fill}} l "
    r"@{\extracolsep{\fill}} l @{\extracolsep{\fill}} l }"
)


@dataclass(frozen=True)
class TableColumn:
    """One logical molecule-table column."""

    header: str
    anchor: str
    molecules: list[Molecule]


def endinput(value: object) -> str:
    """Return a LaTeX input fragment terminated with ``\\endinput``."""
    return f"{value}\\endinput"


def write_fragments(fragments: dict[str, str], output_dir: str | Path = ".") -> None:
    """Write filename-to-content LaTeX fragments to ``output_dir``."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    for filename, content in fragments.items():
        (output_path / filename).write_text(content)


def molecule_link(molecule: Molecule) -> str:
    """Return a hyperlinked mhchem formula for a molecule."""
    return rf"\molref{{{molecule.label}}}{{{molecule.table_formula}}}"


def linked_header(anchor: str, text: str) -> str:
    """Return a hyperlinked table-column header."""
    return rf"\hyperref[{anchor}]{{{text}}}"


def first_detection_sort_keys(
    detections: Iterable[Detection],
) -> dict[str, tuple]:
    """Return first-detection sort keys by molecule label."""
    keys = {}
    for detection in detections:
        label = detection.molecule.label
        key = (
            detection.sortdate,
            detection.molecule.label,
            detection.id,
        )
        if label not in keys or key < keys[label]:
            keys[label] = key
    return keys


def molecules_by_first_detection(
    molecules: Iterable[Molecule],
    detections: Iterable[Detection],
) -> list[Molecule]:
    """Sort molecules by first detection represented in ``detections``."""
    keys = first_detection_sort_keys(detections)
    return sorted(
        molecules,
        key=lambda molecule: keys.get(
            molecule.label,
            (date.max, molecule.label, ""),
        ),
    )


def split_legacy_pair_column(
    values: list[Molecule],
) -> tuple[list[Molecule], list[Molecule]]:
    """Split legacy paired table columns using the historical midpoint rule."""
    midpoint = int(len(values) / 2) + 1
    return values[:midpoint], values[midpoint:]


def split_balanced_columns(
    values: list[Molecule],
    max_rows: int,
) -> list[list[Molecule]]:
    """Split values into balanced chunks with no chunk over ``max_rows``."""
    if not values:
        return [[]]
    nchunks = ceil(len(values) / max_rows)
    chunk_size = ceil(len(values) / nchunks)
    return [
        values[index : index + chunk_size]
        for index in range(0, len(values), chunk_size)
    ]


def table_rows(columns: list[list[Molecule]]) -> list[str]:
    """Render molecule columns as LaTeX table rows."""
    nlines = max(len(column) for column in columns)
    rows = []
    for index in range(nlines):
        row = []
        for column in columns:
            row.append(molecule_link(column[index]) if index < len(column) else "")
        rows.append("\t&\t".join(row) + r"\\")
    return rows


def table_column_rows(columns: list[TableColumn]) -> list[str]:
    """Render logical table columns as LaTeX rows."""
    return table_rows([column.molecules for column in columns])


def table_column_header(columns: list[TableColumn]) -> str:
    """Render logical table-column headers with grouped split columns."""
    cells = []
    index = 0
    while index < len(columns):
        column = columns[index]
        span = 1
        while (
            index + span < len(columns)
            and columns[index + span].header == column.header
            and columns[index + span].anchor == column.anchor
        ):
            span += 1

        linked = linked_header(column.anchor, column.header)
        if span == 1:
            cells.append(linked)
        else:
            cells.append(rf"\multicolumn{{{span}}}{{c}}{{{linked}}}")
        index += span

    return " & ".join(cells) + r" \\"


def tabular_spec(ncols: int) -> str:
    """Return a flexible full-width tabular specification."""
    cols = " ".join(["l"] * ncols)
    return rf"\begin{{tabular*}}{{\textwidth}}{{@{{\extracolsep{{\fill}}}} {cols} @{{}}}}"


def latex_table_fragment(
    *,
    caption: str,
    tabular_spec: str,
    header: str,
    rows: list[str],
    label: str,
) -> str:
    """Render a complete LaTeX table fragment."""
    lines = [
        MOLREF_COMMAND,
        r"\begin{table*}",
        r"\centering",
        rf"\caption{{{caption}}}",
        tabular_spec,
        r"\hline\hline",
        header,
        r"\hline",
        *rows,
        r"\hline",
        r"\end{tabular*}",
        rf"\label{{{label}}}",
        r"\end{table*}\endinput",
    ]
    return "\n".join(lines)


def molecule_count(view: CensusView) -> int:
    """Number of accepted ISM/CSM molecules in the view."""
    return len(view.ism_molecules())


def element_count(view: CensusView) -> int:
    """Number of unique elements in accepted ISM/CSM molecules."""
    elements = set()
    for molecule in view.ism_molecules():
        for element, count in molecule.atom_counts.items():
            if count > 0 and element != "D":
                elements.add(element)
    return len(elements)


def context_molecule_count(
    view: CensusView,
    detection_type: str,
    *,
    include_tentative: bool = False,
    include_disputed: bool = False,
    include_isotopologues: bool = True,
) -> int:
    """Count unique molecules in one detection context."""
    molecules = view.context_molecules(
        detection_type,
        include_tentative=include_tentative,
        include_disputed=include_disputed,
        include_isotopologues=include_isotopologues,
    )
    return len(molecules)


def ppd_molecule_count(view: CensusView) -> int:
    """Number of non-isotopologue molecules detected in PPDs."""
    return context_molecule_count(
        view,
        "ppd",
        include_isotopologues=False,
    )


def ppd_isotopologue_count(view: CensusView) -> int:
    """Number of isotopologue entries detected in PPDs."""
    return len(
        [
            molecule
            for molecule in view.ppd_molecules(include_isotopologues=True)
            if molecule.isotopologue_of is not None
        ]
    )


def exgal_molecule_count(view: CensusView) -> int:
    """Number of accepted extragalactic molecules."""
    return context_molecule_count(view, "exgal")


def exgal_percent(view: CensusView) -> int:
    """Accepted extragalactic molecules as a percentage of ISM/CSM molecules."""
    return round(100 * exgal_molecule_count(view) / molecule_count(view))


def exoplanet_molecule_count(view: CensusView) -> int:
    """Number of accepted exoplanet-atmosphere molecules."""
    return context_molecule_count(view, "exo")


def ice_molecule_count(view: CensusView) -> int:
    """Number of accepted interstellar-ice molecules."""
    return context_molecule_count(view, "ice")


def molecule_detections_by_label(
    detections: Iterable[Detection],
) -> dict[str, list[Detection]]:
    """Group detections by molecule label."""
    grouped: dict[str, list[Detection]] = {}
    for detection in detections:
        grouped.setdefault(detection.molecule.label, []).append(detection)
    return grouped


def radio_percent(view: CensusView) -> int:
    """Percentage of ISM/CSM molecules detected at radio wavelengths."""
    grouped = molecule_detections_by_label(view.ism_detections())
    radio = 0
    for detections in grouped.values():
        if any(
            RADIO_WAVELENGTHS.intersection(detection.wavelengths)
            for detection in detections
        ):
            radio += 1
    return round(100 * radio / molecule_count(view))


def facility_count(view: CensusView) -> int:
    """Number of facilities contributing to accepted ISM/CSM detections."""
    return len(view.facility_counts())


def hydrocarbon_molecules(view: CensusView) -> list[Molecule]:
    """ISM/CSM molecules with C and H and a computable degree of unsaturation."""
    molecules = []
    for molecule in view.ism_molecules():
        atoms = molecule.atom_counts
        if (
            molecule.du is not None
            and atoms.get("C", 0) > 0
            and atoms.get("H", 0) > 0
        ):
            molecules.append(molecule)
    return molecules


def saturated_hydrocarbons(view: CensusView) -> list[Molecule]:
    """Hydrocarbon molecules with zero degree of unsaturation."""
    return [molecule for molecule in hydrocarbon_molecules(view) if molecule.du == 0]


def unsaturated_hydrocarbons(view: CensusView) -> list[Molecule]:
    """Hydrocarbon molecules with nonzero degree of unsaturation."""
    return [molecule for molecule in hydrocarbon_molecules(view) if molecule.du > 0]


def saturated_hydrocarbon_count(view: CensusView) -> int:
    """Number of saturated hydrocarbon molecules."""
    return len(saturated_hydrocarbons(view))


def saturated_hydrocarbon_percent(view: CensusView) -> int:
    """Saturated hydrocarbons as a percentage of hydrocarbon molecules."""
    saturated = saturated_hydrocarbons(view)
    unsaturated = unsaturated_hydrocarbons(view)
    total = len(saturated) + len(unsaturated)
    if total == 0:
        return 0
    return round(100 * len(saturated) / total)


def unsaturated_hydrocarbon_percent(view: CensusView) -> int:
    """Unsaturated hydrocarbons as a percentage of hydrocarbon molecules."""
    saturated = saturated_hydrocarbons(view)
    unsaturated = unsaturated_hydrocarbons(view)
    total = len(saturated) + len(unsaturated)
    if total == 0:
        return 0
    return round(100 * len(unsaturated) / total)


def saturated_hydrocarbon_list(view: CensusView) -> str:
    """LaTeX list of saturated hydrocarbon formulas."""
    molecules = saturated_hydrocarbons(view)
    formulas = [molecule.table_formula for molecule in molecules]
    if not formulas:
        return ""
    if len(formulas) == 1:
        return rf"\ce{{{formulas[0]}}}"
    prefix = ", ".join(rf"\ce{{{formula}}}" for formula in formulas[:-1])
    return prefix + rf", and \ce{{{formulas[-1]}}}"


def molecule_source_types(view: CensusView) -> dict[str, set[str]]:
    """Map ISM/CSM molecule labels to source types represented in detections."""
    source_types: dict[str, set[str]] = {}
    for detection in view.ism_detections():
        values = source_types.setdefault(detection.molecule.label, set())
        values.update(source.type for source in detection.sources)
    return source_types


def radical_percent_by_source_type(view: CensusView, source_type: str) -> int:
    """Percentage of molecules detected in a source type that are radicals."""
    source_types = molecule_source_types(view)
    molecules = {
        molecule.label: molecule
        for molecule in view.ism_molecules()
    }
    in_source = [
        molecules[label]
        for label, types in source_types.items()
        if source_type in types and label in molecules
    ]
    if not in_source:
        return 0
    radicals = [molecule for molecule in in_source if molecule.radical]
    return round(100 * len(radicals) / len(in_source))


def first_detection_years(view: CensusView) -> list[int]:
    """First accepted ISM/CSM detection year for each molecule."""
    grouped = molecule_detections_by_label(view.ism_detections())
    return [
        min(detection.year for detection in detections)
        for detections in grouped.values()
    ]


def detection_rate_since(view: CensusView, start_year: int) -> float:
    """Linear cumulative-detection rate since ``start_year``."""
    years = first_detection_years(view)
    first_year = min(years)
    last_year = max(years)
    x_values = np.arange(first_year, last_year + 1)
    y_values = np.array([
        sum(year <= x_year for year in years)
        for x_year in x_values
    ])
    start_index = np.argwhere(x_values == start_year)[0][0]
    return float(
        np.polynomial.polynomial.Polynomial.fit(
            x_values[start_index:],
            y_values[start_index:],
            1,
        ).convert().coef[1]
    )


def scalar_fragments(view: CensusView) -> dict[str, str]:
    """Return standard scalar LaTeX input fragments for a census view."""
    return {
        "ndetects.tex": endinput(molecule_count(view)),
        "nelems.tex": endinput(element_count(view)),
        "nppds.tex": endinput(ppd_molecule_count(view)),
        "nppdisos.tex": endinput(ppd_isotopologue_count(view)),
        "nexgal.tex": endinput(exgal_molecule_count(view)),
        "nexgalpercent.tex": endinput(exgal_percent(view)),
        "nexos.tex": endinput(exoplanet_molecule_count(view)),
        "nices.tex": endinput(ice_molecule_count(view)),
        "radiopercent.tex": endinput(radio_percent(view)),
        "nscopes.tex": endinput(facility_count(view)),
        "unsatpercent.tex": endinput(unsaturated_hydrocarbon_percent(view)),
        "satlist.tex": endinput(saturated_hydrocarbon_list(view)),
        "nsats.tex": endinput(saturated_hydrocarbon_count(view)),
        "satpercent.tex": endinput(saturated_hydrocarbon_percent(view)),
        "sfr_rad_percent.tex": endinput(
            radical_percent_by_source_type(view, "SFR")
        ),
        "dark_rad_percent.tex": endinput(
            radical_percent_by_source_type(view, "Dark Cloud")
        ),
        "rate_since_1968.tex": endinput(f"{detection_rate_since(view, 1968):.1f}"),
        "rate_since_2005.tex": endinput(f"{detection_rate_since(view, 2005):.1f}"),
    }


def write_scalar_fragments(
    view: CensusView,
    output_dir: str | Path = ".",
) -> dict[str, str]:
    """Write standard scalar LaTeX fragments and return generated content."""
    fragments = scalar_fragments(view)
    write_fragments(fragments, output_dir)
    return fragments


def ism_table_molecules(view: CensusView) -> list[Molecule]:
    """Return non-isotopologue ISM/CSM molecules for the main molecule table."""
    return [
        molecule
        for molecule in view.ism_molecules()
        if molecule.isotopologue_of is None
    ]


def ism_table_detections(view: CensusView) -> list[Detection]:
    """Return non-isotopologue ISM/CSM detections for table ordering."""
    return [
        detection
        for detection in view.ism_detections()
        if detection.molecule.isotopologue_of is None
    ]


def ism_table_columns(
    view: CensusView,
) -> tuple[list[list[Molecule]], list[list[Molecule]]]:
    """Return legacy ISM/CSM table columns for a census view."""
    molecules = molecules_by_first_detection(
        ism_table_molecules(view),
        ism_table_detections(view),
    )
    by_natoms = {
        natoms: [
            molecule
            for molecule in molecules
            if molecule.natoms == natoms
        ]
        for natoms in range(2, 14)
    }

    two_a, two_b = split_legacy_pair_column(by_natoms[2])
    three_a, three_b = split_legacy_pair_column(by_natoms[3])
    four_a, four_b = split_legacy_pair_column(by_natoms[4])
    five_a, five_b = split_legacy_pair_column(by_natoms[5])

    two_seven = [
        two_a,
        two_b,
        three_a,
        three_b,
        four_a,
        four_b,
        five_a,
        five_b,
        by_natoms[6],
        by_natoms[7],
    ]
    eight_more = [
        by_natoms[8],
        by_natoms[9],
        by_natoms[10],
        by_natoms[11],
        by_natoms[12],
        by_natoms[13],
        [molecule for molecule in molecules if molecule.pah],
        [molecule for molecule in molecules if molecule.fullerene],
    ]
    return two_seven, eight_more


def balanced_ism_table_columns(
    view: CensusView,
    *,
    max_rows: int = BALANCED_ISM_MAX_ROWS,
    max_columns: int = BALANCED_ISM_MAX_COLUMNS,
) -> list[list[TableColumn]]:
    """Return balanced ISM/CSM table columns packed into table groups."""
    if max_rows < 1:
        raise ValueError("max_rows must be at least 1.")
    if max_columns < 1:
        raise ValueError("max_columns must be at least 1.")

    molecules = molecules_by_first_detection(
        ism_table_molecules(view),
        ism_table_detections(view),
    )
    categories = []
    for natoms in range(2, 13):
        categories.append(
            (
                f"{natoms} Atoms",
                f"{natoms}atoms",
                [
                    molecule
                    for molecule in molecules
                    if molecule.natoms == natoms
                    and not molecule.pah
                    and not molecule.fullerene
                ],
            )
        )
    categories.append(
        (
            "13+ Atoms",
            "13plusatoms",
            [
                molecule
                for molecule in molecules
                if molecule.natoms >= 13
                and not molecule.pah
                and not molecule.fullerene
            ],
        )
    )
    categories.extend(
        [
            ("PAHs", "pahs", [molecule for molecule in molecules if molecule.pah]),
            (
                "Fullerenes",
                "fullerenes",
                [molecule for molecule in molecules if molecule.fullerene],
            ),
        ]
    )

    category_groups = []
    for header, anchor, values in categories:
        if not values:
            continue
        category_groups.append(
            [
                TableColumn(header, anchor, chunk)
                for chunk in split_balanced_columns(values, max_rows)
            ]
        )

    groups = []
    current_group = []
    for category_group in category_groups:
        remaining = category_group
        while remaining:
            available = max_columns - len(current_group)
            if available == 0:
                groups.append(current_group)
                current_group = []
                available = max_columns

            if len(remaining) <= available:
                current_group.extend(remaining)
                remaining = []
            else:
                if current_group:
                    groups.append(current_group)
                    current_group = []
                else:
                    groups.append(remaining[:max_columns])
                    remaining = remaining[max_columns:]

    if current_group:
        groups.append(current_group)
    return groups


def legacy_ism_table_fragments(
    view: CensusView,
    basename: str = "ism_table",
) -> dict[str, str]:
    """Return the two legacy 2021-style ISM/CSM molecule table fragments."""
    two_seven, eight_more = ism_table_columns(view)
    caption = (
        "List of detected interstellar molecules with {atom_range} atoms, "
        "categorized by number of atoms, and vertically ordered by detection "
        "year.  Column headers and molecule formulas are in-document "
        "hyperlinks in most PDF viewers."
    )
    return {
        f"{basename}_2-7.tex": latex_table_fragment(
            caption=caption.format(atom_range="two to seven"),
            tabular_spec=ISM_TABLE_TWO_SEVEN_SPEC,
            header=(
                r"\multicolumn{2}{c}{\hyperref[2atoms]{2 Atoms}} &"
                r"\multicolumn{2}{c}{\hyperref[3atoms]{3 Atoms}}& "
                r"\multicolumn{2}{c}{\hyperref[4atoms]{4 Atoms}} & "
                r"\multicolumn{2}{c}{\hyperref[5atoms]{5 Atoms}} & "
                r"\hyperref[6atoms]{6 Atoms} & "
                r"\hyperref[7atoms]{7 Atoms} \\"
            ),
            rows=table_rows(two_seven),
            label="two_seven",
        ),
        f"{basename}_8+.tex": latex_table_fragment(
            caption=caption.format(atom_range="eight or more"),
            tabular_spec=ISM_TABLE_EIGHT_MORE_SPEC,
            header=(
                r"\hyperref[8atoms]{8 Atoms} & "
                r"\hyperref[9atoms]{9 Atoms} & "
                r"\hyperref[10atoms]{10 Atoms} & "
                r"\hyperref[11atoms]{11 Atoms} & "
                r"\hyperref[12atoms]{12 Atoms} & "
                r"\hyperref[13atoms]{13 Atoms} & "
                r"\hyperref[pahs]{PAHs} & "
                r"\hyperref[fullerenes]{Fullerenes}  \\"
            ),
            rows=table_rows(eight_more),
            label="eight_more",
        ),
    }


def balanced_ism_table_fragments(
    view: CensusView,
    basename: str = "ism_table",
    *,
    max_rows: int = BALANCED_ISM_MAX_ROWS,
    max_columns: int = BALANCED_ISM_MAX_COLUMNS,
) -> dict[str, str]:
    """Return balanced, future-oriented ISM/CSM molecule table fragments."""
    groups = balanced_ism_table_columns(
        view,
        max_rows=max_rows,
        max_columns=max_columns,
    )
    fragments = {}
    for index, columns in enumerate(groups, start=1):
        fragments[f"{basename}_{index}.tex"] = latex_table_fragment(
            caption=(
                "List of detected interstellar molecules, categorized by "
                "number of atoms or molecular family, and vertically ordered "
                "by detection year.  Column headers and molecule formulas are "
                "in-document hyperlinks in most PDF viewers."
            ),
            tabular_spec=tabular_spec(len(columns)),
            header=table_column_header(columns),
            rows=table_column_rows(columns),
            label=f"ism_molecules_{index}",
        )
    return fragments


def ism_table_fragments(
    view: CensusView,
    basename: str = "ism_table",
    *,
    layout: str = "balanced",
    max_rows: int = BALANCED_ISM_MAX_ROWS,
    max_columns: int = BALANCED_ISM_MAX_COLUMNS,
) -> dict[str, str]:
    """Return ISM/CSM molecule table fragments.

    ``layout="legacy"`` reproduces the two historical audit tables.
    ``layout="balanced"`` produces future-oriented tables constrained by
    ``max_rows`` and ``max_columns``.
    """
    if layout == "legacy":
        return legacy_ism_table_fragments(view, basename)
    if layout == "balanced":
        return balanced_ism_table_fragments(
            view,
            basename,
            max_rows=max_rows,
            max_columns=max_columns,
        )
    raise ValueError("ISM table layout must be 'legacy' or 'balanced'.")


def write_ism_tables(
    view: CensusView,
    output_dir: str | Path = ".",
    basename: str = "ism_table",
    *,
    layout: str = "balanced",
    max_rows: int = BALANCED_ISM_MAX_ROWS,
    max_columns: int = BALANCED_ISM_MAX_COLUMNS,
) -> dict[str, str]:
    """Write standard ISM/CSM molecule tables and return generated content."""
    fragments = ism_table_fragments(
        view,
        basename,
        layout=layout,
        max_rows=max_rows,
        max_columns=max_columns,
    )
    write_fragments(fragments, output_dir)
    return fragments
