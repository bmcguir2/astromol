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
EXGAL_FIRST_TABLE_MAX_PAIRS = 4
EXGAL_TABLE_MAX_PAIRS = 5
PPD_FIRST_TABLE_MAX_PAIRS = 5
PPD_TABLE_MAX_PAIRS = 5
RATE_BY_ATOMS_ONSET_YEARS = {
    2: 1968,
    3: 1968,
    4: 1968,
    5: 1971,
    6: 1970,
    7: 1973,
    8: 1975,
    9: 1974,
    10: 2001,
    11: 2004,
    12: 2001,
    "13+": 2018,
    "PAHs": 2021,
    "Fullerenes": 2010,
}


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


@dataclass(frozen=True)
class DetectionRateFit:
    """Linear detection-rate fit for one molecule-size category."""

    label: str
    onset_year: int
    slope: float
    r_value: float
    r_squared: float


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


def detection_molecule_link(detection: Detection) -> str:
    """Return a linked molecule formula with a tentative marker if needed."""
    link = molecule_link(detection.molecule)
    if detection.status == "tentative":
        return link + r"$^{\dagger}$"
    return link


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


def atom_reference_table_header(atom_counts: Iterable[int]) -> list[str]:
    """Return grouped atom-count and species/reference table headers."""
    atom_counts = list(atom_counts)
    grouped = " & ".join(
        rf"\multicolumn{{2}}{{c}}{{{linked_header(f'{natoms}atoms', f'{natoms} Atoms')}}}"
        for natoms in atom_counts
    )
    labels = " & ".join(["Species & Ref."] * len(atom_counts))
    return [grouped + r" \\", labels + r" \\"]


def tabular_spec(ncols: int) -> str:
    """Return a flexible full-width tabular specification."""
    cols = " ".join(["l"] * ncols)
    return rf"\begin{{tabular*}}{{\textwidth}}{{@{{\extracolsep{{\fill}}}} {cols} @{{}}}}"


def species_reference_tabular_spec(npairs: int) -> str:
    """Return a full-width tabular spec for species/reference column pairs."""
    cols = " ".join(["l l"] * npairs)
    return rf"\begin{{tabular*}}{{\textwidth}}{{@{{\extracolsep{{\fill}}}} {cols} @{{}}}}"


def single_reference_tabular_spec(width: str = r"\columnwidth") -> str:
    """Return a tabular spec for one species/reference column pair."""
    return rf"\begin{{tabular*}}{{{width}}}{{@{{\extracolsep{{\fill}}}} l l @{{}}}}"


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


def rate_fit_end_year(view: CensusView, detections: Iterable[Detection]) -> int:
    """Return the final year to use for rate-table fits."""
    if view.is_current:
        return date.today().year
    return view.census_year or max(detection.year for detection in detections)


def category_detection_years(
    detections: Iterable[Detection],
    category: int | str,
) -> list[int]:
    """Return first-detection years for one atom-count rate category."""
    years = []
    for detection in detections:
        molecule = detection.molecule
        if isinstance(category, int) and molecule.natoms == category:
            years.append(detection.year)
        elif (
            category == "13+"
            and molecule.natoms >= 13
            and not molecule.pah
            and not molecule.fullerene
        ):
            years.append(detection.year)
        elif category == "PAHs" and molecule.pah:
            years.append(detection.year)
        elif category == "Fullerenes" and molecule.fullerene:
            years.append(detection.year)
    return years


def cumulative_counts_by_year(
    detection_years: Iterable[int],
    *,
    onset_year: int,
    end_year: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Return years and cumulative detection counts since the onset year."""
    years = np.arange(onset_year, end_year + 1)
    detection_years = list(detection_years)
    counts = np.array(
        [
            sum(onset_year <= detection_year <= year for detection_year in detection_years)
            for year in years
        ]
    )
    return years, counts


def linear_rate_fit(
    detection_years: Iterable[int],
    *,
    onset_year: int,
    end_year: int,
) -> tuple[float, float, float] | None:
    """Return slope, Pearson R, and R^2 for cumulative detections since onset."""
    detection_years = list(detection_years)
    if not detection_years:
        return None
    years, counts = cumulative_counts_by_year(
        detection_years,
        onset_year=onset_year,
        end_year=end_year,
    )
    if len(years) < 2:
        return None

    slope, intercept = np.polyfit(years, counts, 1)
    if float(np.sum((counts - np.mean(counts)) ** 2)) == 0.0:
        r_value = 0.0
    else:
        r_value = float(np.corrcoef(years, counts)[0, 1])
    return float(slope), r_value, r_value**2


def rate_by_atoms_fits(view: CensusView) -> list[DetectionRateFit]:
    """Return detection-rate fits by atom-count category for ISM/CSM molecules."""
    first_detections = list(
        first_context_detection_by_molecule(view.ism_detections()).values()
    )
    end_year = rate_fit_end_year(view, first_detections)
    fits = []
    for category, onset_year in RATE_BY_ATOMS_ONSET_YEARS.items():
        fit = linear_rate_fit(
            category_detection_years(first_detections, category),
            onset_year=onset_year,
            end_year=end_year,
        )
        if fit is None:
            continue
        slope, r_value, r_squared = fit
        fits.append(
            DetectionRateFit(
                label=str(category),
                onset_year=onset_year,
                slope=slope,
                r_value=r_value,
                r_squared=r_squared,
            )
        )
    return fits


def rate_by_atoms_table_fragment(view: CensusView) -> str:
    """Return the detection-rate-by-atoms table as a LaTeX fragment."""
    lines = [
        r"\begin{table}[htb!]",
        r"\centering",
        (
            r"\caption{Rates (\emph{m}) of detection of new molecules per "
            r"year, sorted by number of atoms per molecule derived from "
            r"linear fits to the data shown in Figure~\ref{cumulative_by_atoms} "
            r"as well as the $R^2$ values of the fits. The start year was "
            r"chosen by the visual onset of a steady detection rate, and is "
            r"given for each fit. Rates and $R^2$ values are obtained from "
            r"least-squares linear fits using NumPy.}"
        ),
        (
            r"\begin{tabular*}{\columnwidth}{c @{\extracolsep{\fill}} "
            r"c @{\extracolsep{\fill}} c @{\extracolsep{\fill}} c }"
        ),
        r"\hline\hline",
        r"\# Atoms    &   \emph{m} (yr$^{-1}$)    &   $R^2$       &   Onset Year      \\",
        r"\hline",
    ]
    for fit in rate_by_atoms_fits(view):
        lines.append(
            f"{fit.label}\t&\t{fit.slope:.2f}\t&\t{fit.r_squared:.2f}"
            f"\t&\t{fit.onset_year}" + r"\\"
        )
    lines.extend(
        [
            r"\hline",
            r"\end{tabular*}",
            r"\label{rates_by_atoms_table}",
            r"\end{table}\endinput",
        ]
    )
    return "\n".join(lines)


def rate_by_atoms_table_fragments(
    view: CensusView,
    filename: str = "rates_by_atoms_table.tex",
) -> dict[str, str]:
    """Return the detection-rate-by-atoms table fragment."""
    return {filename: rate_by_atoms_table_fragment(view)}


def write_rate_by_atoms_table(
    view: CensusView,
    output_dir: str | Path = ".",
    filename: str = "rates_by_atoms_table.tex",
) -> dict[str, str]:
    """Write the detection-rate-by-atoms table and return content."""
    fragments = rate_by_atoms_table_fragments(view, filename)
    write_fragments(fragments, output_dir)
    return fragments


def sorted_count_entries(counts: dict[str, int]) -> list[tuple[str, int]]:
    """Return count entries sorted descending by count, then by label."""
    return sorted(counts.items(), key=lambda item: (-item[1], item[0]))


def paired_count_table_rows(entries: list[tuple[str, int]]) -> list[str]:
    """Render count entries into two side-by-side table-column pairs."""
    nrows = ceil(len(entries) / 2)
    rows = []
    for index in range(nrows):
        left_label, left_count = entries[index]
        if index + nrows < len(entries):
            right_label, right_count = entries[index + nrows]
            rows.append(
                f"{left_label}\t&\t{left_count}\t&\t{right_label}"
                f"\t&\t{right_count}\t" + r"\\"
            )
        else:
            rows.append(
                f"{left_label}\t&\t{left_count}\t&\t\t&\t\t" + r"\\"
            )
    return rows


def paired_count_table_fragment(
    *,
    caption: str,
    label: str,
    item_header: str,
    count_header: str,
    entries: list[tuple[str, int]],
) -> str:
    """Return a two-column-pair LaTeX count table fragment."""
    lines = [
        r"\begin{table}[htb!]",
        r"\centering",
        rf"\caption{{{caption}}}",
        (
            r"\begin{tabular*}{\columnwidth}{l @{\extracolsep{\fill}} "
            r"c @{\extracolsep{\fill}} l @{\extracolsep{\fill}} c }"
        ),
        r"\hline\hline",
        rf"{item_header}	&	{count_header} 	&	{item_header}	&	{count_header} \\",
        r"\hline",
        *paired_count_table_rows(entries),
        r"\hline",
        r"\end{tabular*}",
        rf"\label{{{label}}}",
        r"\end{table}\endinput",
    ]
    return "\n".join(lines)


def facility_table_entries(view: CensusView) -> list[tuple[str, int]]:
    """Return observing-facility contribution counts for the facility table."""
    return sorted_count_entries(view.facility_counts(key="latex_name"))


def facility_table_fragment(view: CensusView) -> str:
    """Return the observing-facility count table as a LaTeX fragment."""
    return paired_count_table_fragment(
        caption=r"Total number of detections for each facility listed in \S\ref{known}.",
        label="detects_by_scope",
        item_header="Facility",
        count_header=r"\#",
        entries=facility_table_entries(view),
    )


def facility_table_fragments(
    view: CensusView,
    filename: str = "facilities_table.tex",
) -> dict[str, str]:
    """Return the observing-facility count table fragment."""
    return {filename: facility_table_fragment(view)}


def write_facility_table(
    view: CensusView,
    output_dir: str | Path = ".",
    filename: str = "facilities_table.tex",
) -> dict[str, str]:
    """Write the observing-facility count table and return content."""
    fragments = facility_table_fragments(view, filename)
    write_fragments(fragments, output_dir)
    return fragments


def source_table_entries(view: CensusView) -> list[tuple[str, int]]:
    """Return source contribution counts for the source table."""
    return sorted_count_entries(
        view.source_counts(
            key="latex_name",
            group_diffuse_cloud=True,
            diffuse_cloud_label="Diffuse Cloud",
        )
    )


def source_table_fragment(view: CensusView) -> str:
    """Return the source contribution count table as a LaTeX fragment."""
    return paired_count_table_fragment(
        caption=(
            r"Total number of detections that each source contributed to for "
            r"the molecules listed in \S\ref{known}. Detections made in "
            r"diffuse clouds along the line of sight to a background source "
            r"have been consolidated into `Diffuse Cloud,' and detections in "
            r"closely located regions have been grouped together as well "
            r"(e.g. Sgr B2(OH), Sgr B2(N), Sgr B2(S), and Sgr B2(M) are all "
            r"considered Sgr B2)."
        ),
        label="detects_by_source",
        item_header="Source",
        count_header=r"\#",
        entries=source_table_entries(view),
    )


def source_table_fragments(
    view: CensusView,
    filename: str = "source_table.tex",
) -> dict[str, str]:
    """Return the source contribution count table fragment."""
    return {filename: source_table_fragment(view)}


def write_source_table(
    view: CensusView,
    output_dir: str | Path = ".",
    filename: str = "source_table.tex",
) -> dict[str, str]:
    """Write the source contribution count table and return content."""
    fragments = source_table_fragments(view, filename)
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


def first_context_detection_by_molecule(
    detections: Iterable[Detection],
) -> dict[str, Detection]:
    """Return the first secure detection per molecule, falling back to tentative."""
    grouped = molecule_detections_by_label(detections)
    selected = {}
    for label, molecule_detections in grouped.items():
        secure = [
            detection
            for detection in molecule_detections
            if detection.status == "secure"
        ]
        candidates = secure or molecule_detections
        selected[label] = min(
            candidates,
            key=lambda detection: (
                detection.sortdate,
                detection.id,
            ),
        )
    return selected


def exgal_table_detections(view: CensusView) -> list[Detection]:
    """Return non-isotopologue exgal detections for the molecule table.

    Secure detections are included by accepted census membership. Tentative
    detections are included by introduced census membership and rendered with a
    dagger marker.
    """
    detections_by_label = first_context_detection_by_molecule(
        view.exgal_detections(include_tentative=True)
    )
    return [
        detections_by_label[molecule.label]
        for molecule in view.db.molecules.values()
        if molecule.label in detections_by_label
        and molecule.isotopologue_of is None
    ]


def chunked_atom_groups(
    atom_counts: list[int],
    *,
    first_group_size: int,
    group_size: int,
) -> list[tuple[int, ...]]:
    """Split atom-count labels into first-page and continuation groups."""
    if not atom_counts:
        return []

    groups = [tuple(atom_counts[:first_group_size])]
    remaining = atom_counts[first_group_size:]
    groups.extend(
        tuple(remaining[index : index + group_size])
        for index in range(0, len(remaining), group_size)
    )
    return [group for group in groups if group]


def exgal_table_atom_groups(view: CensusView) -> list[tuple[int, ...]]:
    """Return non-empty atom-count groups for the exgal molecule table."""
    atom_counts = sorted(
        {
            detection.molecule.natoms
            for detection in exgal_table_detections(view)
        }
    )
    return chunked_atom_groups(
        atom_counts,
        first_group_size=EXGAL_FIRST_TABLE_MAX_PAIRS,
        group_size=EXGAL_TABLE_MAX_PAIRS,
    )


def exgal_table_columns(
    view: CensusView,
    atom_groups: Iterable[Iterable[int]] | None = None,
) -> list[list[list[Detection]]]:
    """Return exgal table columns split by configured atom-count groups."""
    detections = exgal_table_detections(view)
    if atom_groups is None:
        atom_groups = exgal_table_atom_groups(view)
    return [
        [
            [
                detection
                for detection in detections
                if detection.molecule.natoms == natoms
            ]
            for natoms in atom_group
        ]
        for atom_group in atom_groups
    ]


def observation_ref_bibcodes(detection: Detection) -> list[str]:
    """Return observation-reference BibTeX keys for a detection."""
    return [ref.bibcode for ref in detection.refs.get("observation", [])]


def reference_numbers(
    bibcodes: Iterable[str],
    references: dict[str, int],
) -> str:
    """Assign and render stable numeric reference labels for BibTeX keys."""
    values = []
    for bibcode in bibcodes:
        if bibcode not in references:
            references[bibcode] = len(references) + 1
        values.append(str(references[bibcode]))
    return ", ".join(values)


def detection_reference_table_rows(
    columns: list[list[Detection]],
    references: dict[str, int],
) -> list[str]:
    """Render rows for molecule/reference table columns."""
    nlines = max(len(column) for column in columns)
    rows = []
    for index in range(nlines):
        row = []
        for column in columns:
            if index < len(column):
                detection = column[index]
                row.extend(
                    [
                        detection_molecule_link(detection),
                        reference_numbers(
                            observation_ref_bibcodes(detection),
                            references,
                        ),
                    ]
                )
            else:
                row.extend(["", ""])
        rows.append("\t&\t".join(row) + r"\\")
    return rows


def single_detection_reference_table_rows(
    detections: Iterable[Detection],
    references: dict[str, int],
) -> list[str]:
    """Render rows for a one-column species/reference table."""
    rows = []
    for detection in detections:
        rows.append(
            "\t&\t".join(
                [
                    detection_molecule_link(detection),
                    reference_numbers(
                        observation_ref_bibcodes(detection),
                        references,
                    ),
                ]
            )
            + r"\\"
        )
    return rows


def numbered_reference_notes(references: dict[str, int]) -> str:
    """Render numbered citation notes from BibTeX key to index mapping."""
    entries = sorted(references.items(), key=lambda item: item[1])
    refs = " ".join(
        rf"[{index}] \citet{{{bibcode}}}"
        for bibcode, index in entries
    )
    return rf"\textbf{{References:}} {refs}\\"


def exgal_table_fragment(view: CensusView) -> str:
    """Return the external-galaxy molecule table as a LaTeX fragment."""
    references: dict[str, int] = {}
    atom_groups = exgal_table_atom_groups(view)
    table_groups = exgal_table_columns(view, atom_groups)
    lines = [
        MOLREF_COMMAND,
        r"\begin{table*}",
        r"\centering",
        (
            r"\caption{List of molecules detected in external galaxies with "
            r"references to the first detections.  Tentative detections are "
            r"indicated, and some extra references are occasionally provided "
            r"for context.}"
        ),
    ]

    for group_index, (atom_counts, columns) in enumerate(
        zip(atom_groups, table_groups)
    ):
        if group_index > 0:
            lines.extend([r"\hline\hline", r"\end{tabular*}"])
        lines.append(species_reference_tabular_spec(len(atom_counts)))
        if group_index == 0:
            lines.append(r"\hline\hline")
        lines.extend(atom_reference_table_header(atom_counts))
        lines.append(r"\hline")
        lines.extend(detection_reference_table_rows(columns, references))

    lines.extend(
        [
            r"\hline",
            r"\end{tabular*}",
            r"\justify",
            r"$^{\dagger}$Tentative detection\\",
            numbered_reference_notes(references),
            r"\label{exgal_mols}",
            r"\end{table*}\endinput",
        ]
    )
    return "\n".join(lines)


def exgal_table_fragments(
    view: CensusView,
    filename: str = "exgal_table.tex",
) -> dict[str, str]:
    """Return the external-galaxy molecule table fragment."""
    return {filename: exgal_table_fragment(view)}


def write_exgal_table(
    view: CensusView,
    output_dir: str | Path = ".",
    filename: str = "exgal_table.tex",
) -> dict[str, str]:
    """Write the external-galaxy molecule table and return generated content."""
    fragments = exgal_table_fragments(view, filename)
    write_fragments(fragments, output_dir)
    return fragments


def ppd_table_detections(view: CensusView) -> list[Detection]:
    """Return PPD detections for the molecule table, including isotopologues.

    The PPD table intentionally includes isotope records. Ordering is
    parent-first: each parent molecule is followed by any detected PPD
    isotopologues that point to it.
    """
    detections_by_label = first_context_detection_by_molecule(
        view.ppd_detections(include_isotopologues=True)
    )
    isotopologues_by_parent: dict[str, list[Detection]] = {}
    for molecule in view.db.molecules.values():
        if (
            molecule.isotopologue_of is not None
            and molecule.label in detections_by_label
        ):
            isotopologues_by_parent.setdefault(
                molecule.isotopologue_of,
                [],
            ).append(detections_by_label[molecule.label])

    ordered = []
    added = set()
    for molecule in view.db.molecules.values():
        if molecule.isotopologue_of is not None:
            continue
        if molecule.label in detections_by_label:
            ordered.append(detections_by_label[molecule.label])
            added.add(molecule.label)
        for detection in isotopologues_by_parent.get(molecule.label, []):
            ordered.append(detection)
            added.add(detection.molecule.label)

    # Keep orphaned isotope records visible instead of silently dropping them.
    for molecule in view.db.molecules.values():
        if molecule.label in detections_by_label and molecule.label not in added:
            ordered.append(detections_by_label[molecule.label])

    return ordered


def ppd_table_atom_groups(view: CensusView) -> list[tuple[int, ...]]:
    """Return non-empty atom-count groups for the PPD molecule table."""
    atom_counts = sorted(
        {
            detection.molecule.natoms
            for detection in ppd_table_detections(view)
        }
    )
    return chunked_atom_groups(
        atom_counts,
        first_group_size=PPD_FIRST_TABLE_MAX_PAIRS,
        group_size=PPD_TABLE_MAX_PAIRS,
    )


def ppd_table_columns(
    view: CensusView,
    atom_groups: Iterable[Iterable[int]] | None = None,
) -> list[list[list[Detection]]]:
    """Return PPD table columns split by configured atom-count groups."""
    detections = ppd_table_detections(view)
    if atom_groups is None:
        atom_groups = ppd_table_atom_groups(view)
    return [
        [
            [
                detection
                for detection in detections
                if detection.molecule.natoms == natoms
            ]
            for natoms in atom_group
        ]
        for atom_group in atom_groups
    ]


def ppd_table_fragment(view: CensusView) -> str:
    """Return the protoplanetary-disk molecule table as a LaTeX fragment."""
    references: dict[str, int] = {}
    atom_groups = ppd_table_atom_groups(view)
    table_groups = ppd_table_columns(view, atom_groups)
    lines = [
        MOLREF_COMMAND,
        r"\begin{table*}",
        r"\centering",
        (
            r"\caption{List of molecules, including rare isotopic species, "
            r"detected in protoplanetary disks, with references to "
            r"representative detections.  The earliest reported detection of a "
            r"species in the literature is provided on a best-effort basis.  "
            r"Tentative and disputed detections are not included (see text).}"
        ),
    ]

    for group_index, (atom_counts, columns) in enumerate(
        zip(atom_groups, table_groups)
    ):
        if group_index > 0:
            lines.extend([r"\hline\hline", r"\end{tabular*}"])
        lines.append(species_reference_tabular_spec(len(atom_counts)))
        if group_index == 0:
            lines.append(r"\hline\hline")
        lines.extend(atom_reference_table_header(atom_counts))
        lines.append(r"\hline")
        lines.extend(detection_reference_table_rows(columns, references))

    lines.extend(
        [
            r"\hline",
            r"\end{tabular*}",
            r"\justify",
            numbered_reference_notes(references),
            r"\label{ppd_mols}",
            r"\end{table*}\endinput",
        ]
    )
    return "\n".join(lines)


def ppd_table_fragments(
    view: CensusView,
    filename: str = "ppd_table.tex",
) -> dict[str, str]:
    """Return the protoplanetary-disk molecule table fragment."""
    return {filename: ppd_table_fragment(view)}


def write_ppd_table(
    view: CensusView,
    output_dir: str | Path = ".",
    filename: str = "ppd_table.tex",
) -> dict[str, str]:
    """Write the protoplanetary-disk molecule table and return content."""
    fragments = ppd_table_fragments(view, filename)
    write_fragments(fragments, output_dir)
    return fragments


def parent_first_context_detections(
    db,
    detections: Iterable[Detection],
    *,
    include_isotopologues: bool,
) -> list[Detection]:
    """Order context detections in database order, parent before isotopologues."""
    detections_by_label = first_context_detection_by_molecule(detections)
    isotopologues_by_parent: dict[str, list[Detection]] = {}
    for molecule in db.molecules.values():
        if (
            molecule.isotopologue_of is not None
            and molecule.label in detections_by_label
        ):
            isotopologues_by_parent.setdefault(
                molecule.isotopologue_of,
                [],
            ).append(detections_by_label[molecule.label])

    ordered = []
    added = set()
    for molecule in db.molecules.values():
        if molecule.isotopologue_of is not None:
            continue
        if molecule.label in detections_by_label:
            ordered.append(detections_by_label[molecule.label])
            added.add(molecule.label)
        if include_isotopologues:
            for detection in isotopologues_by_parent.get(molecule.label, []):
                ordered.append(detection)
                added.add(detection.molecule.label)

    if include_isotopologues:
        for molecule in db.molecules.values():
            if molecule.label in detections_by_label and molecule.label not in added:
                ordered.append(detections_by_label[molecule.label])

    return ordered


def exoplanet_table_detections(
    view: CensusView,
    *,
    include_isotopologues: bool = False,
) -> list[Detection]:
    """Return exoplanet detections for the molecule table."""
    return parent_first_context_detections(
        view.db,
        view.exoplanet_detections(include_isotopologues=include_isotopologues),
        include_isotopologues=include_isotopologues,
    )


def exoplanet_table_fragment(
    view: CensusView,
    *,
    include_isotopologues: bool = False,
) -> str:
    """Return the exoplanet-atmosphere molecule table as a LaTeX fragment."""
    references: dict[str, int] = {}
    detections = exoplanet_table_detections(
        view,
        include_isotopologues=include_isotopologues,
    )
    lines = [
        MOLREF_COMMAND,
        r"\begin{table}",
        r"\centering",
        (
            r"\caption{List of molecules detected in exoplanetary atmospheres, "
            r"with references to representative detections.  Tentative and "
            r"disputed detections are not included.}"
        ),
        single_reference_tabular_spec(),
        r"\hline\hline",
        r"Species & References\\",
        r"\hline",
        *single_detection_reference_table_rows(detections, references),
        r"\hline",
        r"\end{tabular*}",
        r"\justify",
        numbered_reference_notes(references),
        r"\label{exoplanet_mols}",
        r"\end{table}\endinput",
    ]
    return "\n".join(lines)


def exoplanet_table_fragments(
    view: CensusView,
    filename: str = "exo_table.tex",
    *,
    include_isotopologues: bool = False,
) -> dict[str, str]:
    """Return the exoplanet-atmosphere molecule table fragment."""
    return {
        filename: exoplanet_table_fragment(
            view,
            include_isotopologues=include_isotopologues,
        )
    }


def write_exoplanet_table(
    view: CensusView,
    output_dir: str | Path = ".",
    filename: str = "exo_table.tex",
    *,
    include_isotopologues: bool = False,
) -> dict[str, str]:
    """Write the exoplanet-atmosphere table and return content."""
    fragments = exoplanet_table_fragments(
        view,
        filename,
        include_isotopologues=include_isotopologues,
    )
    write_fragments(fragments, output_dir)
    return fragments


def ice_table_detections(
    view: CensusView,
    *,
    include_tentative: bool = True,
    include_isotopologues: bool = False,
) -> list[Detection]:
    """Return ice detections for the molecule table."""
    return parent_first_context_detections(
        view.db,
        view.ice_detections(
            include_tentative=include_tentative,
            include_isotopologues=include_isotopologues,
        ),
        include_isotopologues=include_isotopologues,
    )


def ice_table_fragment(
    view: CensusView,
    *,
    include_tentative: bool = True,
    include_isotopologues: bool = False,
) -> str:
    """Return the interstellar-ice molecule table as a LaTeX fragment."""
    references: dict[str, int] = {}
    detections = ice_table_detections(
        view,
        include_tentative=include_tentative,
        include_isotopologues=include_isotopologues,
    )
    tentative_note = (
        [r"$^{\dagger}$Tentative detection\\"]
        if any(detection.status == "tentative" for detection in detections)
        else []
    )
    lines = [
        MOLREF_COMMAND,
        r"\begin{table}",
        r"\centering",
        (
            r"\caption{List of molecules detected in interstellar ices, with "
            r"references to representative detections. Tentative detections "
            r"are indicated.}"
        ),
        single_reference_tabular_spec(),
        r"\hline\hline",
        r"Species & References\\",
        r"\hline",
        *single_detection_reference_table_rows(detections, references),
        r"\hline",
            r"\end{tabular*}",
            r"\justify",
            *tentative_note,
            numbered_reference_notes(references),
            r"\label{ice_mols}",
            r"\end{table}\endinput",
    ]
    return "\n".join(lines)


def ice_table_fragments(
    view: CensusView,
    filename: str = "ice_table.tex",
    *,
    include_tentative: bool = True,
    include_isotopologues: bool = False,
) -> dict[str, str]:
    """Return the interstellar-ice molecule table fragment."""
    return {
        filename: ice_table_fragment(
            view,
            include_tentative=include_tentative,
            include_isotopologues=include_isotopologues,
        )
    }


def write_ice_table(
    view: CensusView,
    output_dir: str | Path = ".",
    filename: str = "ice_table.tex",
    *,
    include_tentative: bool = True,
    include_isotopologues: bool = False,
) -> dict[str, str]:
    """Write the interstellar-ice table and return content."""
    fragments = ice_table_fragments(
        view,
        filename,
        include_tentative=include_tentative,
        include_isotopologues=include_isotopologues,
    )
    write_fragments(fragments, output_dir)
    return fragments
