"""LaTeX output helpers for census manuscripts."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

import numpy as np

from .census import CensusView
from .models import Detection, Molecule


RADIO_WAVELENGTHS = {"cm", "mm", "sub-mm"}


def endinput(value: object) -> str:
    """Return a LaTeX input fragment terminated with ``\\endinput``."""
    return f"{value}\\endinput"


def write_fragments(fragments: dict[str, str], output_dir: str | Path = ".") -> None:
    """Write filename-to-content LaTeX fragments to ``output_dir``."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    for filename, content in fragments.items():
        (output_path / filename).write_text(content)


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
    )
    if not include_isotopologues:
        molecules = [
            molecule
            for molecule in molecules
            if molecule.isotopologue_of is None
        ]
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
            for molecule in view.ppd_molecules()
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
