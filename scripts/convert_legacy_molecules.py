"""Convert legacy astromol molecule definitions into preview JSON files.

This script intentionally does not overwrite the active data files. It parses
``astromol/data/molecules_legacy.py`` with ``ast`` rather than importing it,
converts only mechanically safe fields, and writes an issue report for data
that requires curation.
"""

from __future__ import annotations

import argparse
import ast
import json
import re
import unicodedata
from collections import Counter, defaultdict
from pathlib import Path

from molmass import Formula


REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = REPO_ROOT / "astromol" / "data"
LEGACY_PATH = DATA_DIR / "molecules_legacy.py"
REFERENCES_PATH = DATA_DIR / "references.bib"
SOURCES_PATH = DATA_DIR / "sources.json"
TELESCOPES_PATH = DATA_DIR / "telescopes.json"

DEFAULT_MOLECULES_OUT = DATA_DIR / "molecules.preview.json"
DEFAULT_DETECTIONS_OUT = DATA_DIR / "detections.preview.json"
DEFAULT_ISSUES_OUT = DATA_DIR / "legacy_conversion_issues.json"
DEFAULT_TRIAGE_OUT = DATA_DIR / "legacy_conversion_triage.md"

REFERENCE_KEY_ALIASES = {
    "Brown:1982ur": "Brown:1982:1747",
    "Brown:1986lp": "Brown:1986:1296",
    "Crabtree:2016fj": "Crabtree:2016:124201",
    "Heineking:1994op": "Heineking:1994:1177",
    "Kaushik:1982ld": "Kaushik:1982:117",
    "Sakaizumi:1976uu": "Sakaizumi:1976:2908",
}

FREE_TEXT_REFERENCE_ALIASES = {
    # Confirmed legacy citation typos or omissions. Keep these in code because
    # references.bib is a Zotero export and should not be hand-edited.
    "Benson & Flygare 1970 J Am Chem Soc 92, 7523": "Flygare:1970:7523",
    "Bouche et al. 1973 J Mol Struct 18, 211": "Bouchy:1973:211",
    "Cernicharo et al. 2021 A&AL 649, 15": "Cernicharo:2021:L15a",
    "Cord et al. 1968 Microwave Spectral Tables V5": "Cord:1968:",
    "Dixon 1959 Can J. Phys. 37, 1171 and Klaus et al. 1997 A&A 322, L1": [
        "Dixon:1959:1171",
        "Klaus:1997:L1",
    ],
    "Fayolle et al. 2017 Nature Astron. 1, 702": "Fayolle:2017:703",
    "Goldhaber and Betz 1977 ApJ 279, L55": "Goldhaber:1984:L55",
    "Gottlieb et al. 200 JCP 113, 1910": "Gottlieb:2000:1910",
    "Haas etal. 1994 JMS 167, 176": "Haas:1994:176",
    "Jevons 1932 Phys Soc. pp 177-179": "Jevons:1932:177",
    "Kasai & Myers JCP 30, 1096": "Kasai:1959:1096",
    "Kattija-Ari & Harmony et al. 1980 International Journal of Quantum Chemistry 18, 443": "Kattija-Ari:2009:443",
    "Kaushik 1977 CPL 49, 90": "Kaushik:1977:89",
    "Kessler et al. Phys Rev 79, 54": "Kessler:1950:54",
    "Kruger et al. 2010 Ang. Chem. 23, 1644": "Guarnieri:1992:39",
    "Lee & Amano 1987 ApJ 323": "Lee:1987:L145",
    "Loomis et al. 2013 ApJL 765, L10": "Loomis:2013:L9",
    "McCarthy et al. JCP 110, 1064": "McCarthy:1999:10645",
    "Miller et al. 1962 JMS 8, 153": "Millen:1962:153",
    "Nakimi et al. 1998 JMS 191, 176": "Namiki:1998:176",
    "Oshima & Endo 1993 JMS 159, 458": "Ohshima:1993:458",
    "Shinegari 1967 J Phys Soc Jpn 23, 404": "Shigenari:1967:404",
    "Steenbeckeliers 1968 Ann. Soc. Sci. Brux 82, 331": "Steenbeckeliers:1968:331",
    "Takano et al. 1986 JMS 26, 157": "Takano:1968:157",
    "Thaddues & Turner 1975 ApJ 201, L25": "Thaddeus:1975:L25",
    "Thomas & Dalby 1968 Can. J. Phys. 46, 2815": "Thomson:1968:2815",
    "Additional work used in Belloche et al. 2019 A&A 628, A10 to be reported in Medvedev et al. in prep as of 9/16/2019.": "Tyree:2022:111706",
    "Zaleski et al. 2013 ApJ 765, L9": "Zaleski:2013:L10",
    "Zukerman et al. 1971 ApJ 163, L41": "Zuckerman:1971:L41",
    "Zukerman et al. 1975 ApJ 196, L99": "Zuckerman:1975:L99",
}

SOURCE_ALIASES = {
    "LOSCloud": "DiffuseCloud",
    "Orion": "OrionKL",
    "TYC 8998-760-1 b": "TYC89987601b",
}

FORMULA_OVERRIDES = {
    "iC4H8": "C4H8",
    "H2^{13}CO": "H2[13C]O",
}

ISOTOPOLOGUE_FORMULA_OVERRIDES = {
    "C15N": "C[15N]",
    "13CO": "[13C]O",
    "C18O": "C[18O]",
    "C17O": "C[17O]",
    "13C17O": "[13C][17O]",
    "13C18O": "[13C][18O]",
    "HD": "H[2H]",
    "C34S": "C[34S]",
    "13CS": "[13C]S",
    "34SO": "[34S]O",
    "33SO": "[33S]O",
    "DCO+": "[2H]CO+",
    "H13CO+": "H[13C]O+",
    "HC18O+": "HC[18O]+",
    "DCN": "[2H]CN",
    "H13CN": "H[13C]N",
    "HC15N": "HC[15N]",
    "DNC": "[2H]NC",
    "N2D+": "N2[2H]+",
    "C2D": "C2[2H]",
    "13CO2": "[13C]O2",
    "H2^{13}CO": "H2[13C]O",
    "H13CCH": "H[13C]CH",
    "CH3D": "CH3[2H]",
}

NAME_OVERRIDES = {
    "H3p": "Trihydrogen ion",
    "HNCCC": "Iminopropadienylidene",
    "HC4N": "cyanoethynylmethylene",
    "HC7NHp": "Protonated-2,4,6-heptatriynenitrile",
}

ALLOWED_WAVELENGTHS = {"cm", "mm", "sub-mm", "IR", "UV", "Vis"}
EXTRA_CONTEXTS = ("ice", "ppd", "exgal", "exo")
NESTED_ISOTOPOLOGUE_FIELDS = ("ice_isos", "ppd_isos", "exgal_isos", "exo_isos")


class LegacyConverter:
    def __init__(self, only: set[str] | None = None):
        self.only = only
        (
            self.reference_ids,
            self.reference_indexes,
            self.reference_years,
        ) = self._load_references()
        self.reference_key_aliases = self._normalized_aliases(REFERENCE_KEY_ALIASES)
        self.free_text_reference_aliases = self._normalized_aliases(
            FREE_TEXT_REFERENCE_ALIASES
        )
        self.source_nicks = self._load_nicks(SOURCES_PATH)
        self.telescope_nicks = self._load_nicks(TELESCOPES_PATH)
        self.reference_mappings = []
        self.issues = []

    def convert(self):
        legacy_molecules = self._load_legacy_molecules()
        molecules = []
        detections = []
        molecule_labels = set()

        for legacy in legacy_molecules:
            if self.only and not self._matches_only(legacy):
                continue

            molecule = self._convert_molecule(legacy)
            molecules.append(molecule)
            molecule_labels.add(molecule["label"])

            for isotopologue in self._convert_nested_isotopologue_molecules(
                legacy, molecule
            ):
                if isotopologue["label"] in molecule_labels:
                    continue
                molecules.append(isotopologue)
                molecule_labels.add(isotopologue["label"])

            detection = self._convert_main_detection(legacy, molecule["label"])
            if detection is not None:
                detections.append(detection)

            detections.extend(self._convert_nested_isotopologue_detections(legacy))

            self._report_omitted_extra_contexts(legacy)

        self._validate_preview(molecules, detections)
        return molecules, detections, self._issue_report(molecules, detections)

    def _load_references(self):
        try:
            import bibtexparser
        except ImportError as exc:
            raise ImportError(
                "Legacy conversion requires the 'bibtexparser' package."
            ) from exc

        with open(REFERENCES_PATH) as handle:
            database = bibtexparser.load(handle)

        reference_ids = {entry["ID"] for entry in database.entries if "ID" in entry}
        reference_years = {}
        indexes = {
            "ads": defaultdict(list),
            "ads_volume": defaultdict(list),
            "text": defaultdict(list),
            "text_volume": defaultdict(list),
            "zotero_base": defaultdict(list),
        }

        for entry in database.entries:
            ref_id = entry.get("ID")
            if not ref_id:
                continue

            year = self._normalize_year(entry.get("year"))
            if year is not None:
                reference_years[ref_id] = int(year)

            volume = self._normalize_volume(entry.get("volume"))
            pages = self._page_variants(entry.get("pages") or entry.get("page"))
            surnames = self._reference_surname_keys(entry)
            initials = {surname[:1].upper() for surname in surnames if surname}

            for page in pages:
                for surname in surnames:
                    indexes["text"][(year, surname, page)].append(ref_id)
                    if volume is not None:
                        indexes["text_volume"][
                            (year, surname, volume, page)
                        ].append(ref_id)
                for initial in initials:
                    indexes["ads"][(year, initial, page)].append(ref_id)
                    if volume is not None:
                        indexes["ads_volume"][
                            (year, initial, volume, page)
                        ].append(ref_id)

            base = self._zotero_base_key(ref_id)
            if base != ref_id:
                indexes["zotero_base"][base].append(ref_id)

        return reference_ids, indexes, reference_years

    def _normalize_year(self, value):
        match = re.search(r"\d{4}", str(value or ""))
        return match.group(0) if match else None

    def _normalize_volume(self, value):
        text = re.sub(r"[^A-Za-z0-9]", "", str(value or "")).upper()
        return text or None

    def _reference_surname_keys(self, entry):
        surnames = set()

        ref_id = entry.get("ID")
        if ref_id and ":" in ref_id:
            surnames.add(self._normalize_lookup_text(ref_id.split(":", 1)[0]))

        authors = entry.get("author") or entry.get("editor")
        if authors:
            surname = self._first_author_surname(authors)
            if surname:
                surnames.add(self._normalize_lookup_text(surname))

        return {surname for surname in surnames if surname}

    def _first_author_surname(self, authors):
        first = authors.split(" and ", 1)[0].strip().strip("{}")
        if "," in first:
            return first.split(",", 1)[0].strip("{} ")
        return first.split()[-1].strip("{} ") if first.split() else None

    def _normalize_lookup_text(self, value):
        text = str(value or "")
        text = re.sub(r"\\['\"`^~=.]\\{?([A-Za-z])\\}?", r"\1", text)
        text = text.replace("{", "").replace("}", "")
        text = unicodedata.normalize("NFKD", text)
        text = "".join(ch for ch in text if not unicodedata.combining(ch))
        return re.sub(r"[^A-Za-z0-9]", "", text).lower()

    def _zotero_base_key(self, ref_id):
        return re.sub(r"(?<=\d|[A-Z])([a-z])$", "", ref_id)

    def _page_variants(self, page):
        if page is None:
            return set()

        first = str(page).split("--", 1)[0].split("-", 1)[0].strip("{} ")
        cleaned = re.sub(r"[^A-Za-z0-9]", "", first).upper()
        return {cleaned} if cleaned else set()

    def _normalized_aliases(self, aliases):
        return {
            self._reference_alias_key(original): mapped
            for original, mapped in aliases.items()
        }

    def _reference_alias_key(self, value):
        return re.sub(r"\s+", " ", str(value or "").strip()).lower()

    def _load_nicks(self, path):
        with open(path) as handle:
            return {entry["nick"] for entry in json.load(handle)}

    def _load_legacy_molecules(self):
        tree = ast.parse(LEGACY_PATH.read_text())
        molecules = []

        for node in tree.body:
            if not self._is_top_level_molecule_assignment(node):
                continue

            var_name = next(
                target.id for target in node.targets if isinstance(target, ast.Name)
            )
            molecule = self._parse_molecule_call(node.value)
            molecule["__legacy_var"] = var_name
            molecule["__line"] = node.lineno
            molecules.append(molecule)

        return molecules

    def _is_top_level_molecule_assignment(self, node):
        return (
            isinstance(node, ast.Assign)
            and isinstance(node.value, ast.Call)
            and isinstance(node.value.func, ast.Name)
            and node.value.func.id == "Molecule"
            and any(isinstance(target, ast.Name) for target in node.targets)
        )

    def _parse_molecule_call(self, call):
        data = {"__line": getattr(call, "lineno", None)}
        for keyword in call.keywords:
            if keyword.arg is not None:
                data[keyword.arg] = self._parse_value(keyword.value)
        return data

    def _parse_value(self, node):
        if isinstance(node, ast.Constant):
            return node.value
        if isinstance(node, ast.Name):
            return node.id
        if isinstance(node, (ast.List, ast.Tuple, ast.Set)):
            return [self._parse_value(value) for value in node.elts]
        if isinstance(node, ast.Dict):
            return {
                self._parse_value(key): self._parse_value(value)
                for key, value in zip(node.keys, node.values)
            }
        if (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Name)
            and node.func.id == "Molecule"
        ):
            return self._parse_molecule_call(node)

        raise ValueError(
            f"Unsupported AST node at line {getattr(node, 'lineno', 'unknown')}: "
            f"{ast.dump(node, include_attributes=False)}"
        )

    def _matches_only(self, legacy):
        candidates = {
            legacy.get("__legacy_var"),
            legacy.get("label"),
            legacy.get("astromol_name"),
            legacy.get("formula"),
            self._molecule_label(legacy),
        }
        return bool(candidates & self.only)

    def _convert_molecule(self, legacy):
        label = self._molecule_label(legacy)
        original_formula = legacy.get("formula")
        formula = self._normalize_formula(original_formula, legacy, label)
        name = self._molecule_name(legacy, formula, label)
        lab_refs = self._refs_from_fields(
            legacy, "l_ref_bib_ids", "l_refs", "lab"
        )

        molecule = {
            "name": name,
            "formula": formula,
            "table_formula": self._table_formula(legacy, original_formula, formula),
            "label": label,
            "note": self._none_if_blank(legacy.get("notes")),
            "iupac_name": None,
            "selfies": None,
            "synonyms": [],
            "smiles": self._none_if_blank(legacy.get("smiles")),
            "canonical_smiles": None,
            "inchi": None,
            "inchikey": None,
            "radical_override": None,
            "fullerene": bool(legacy.get("fullerene", False)),
            "pah": bool(legacy.get("pah", False)),
            "n_rings": 0,
            "cyclic": bool(legacy.get("cyclic", False)),
            "rotcon": self._rotcon(legacy, lab_refs),
            "dipole": self._dipole(legacy, lab_refs),
            "refs": {
                "lab": lab_refs,
                "computation": [],
            },
            "isotopologue_of": None,
            "latex_header": None,
            "latex_notes": None,
        }

        return molecule

    def _convert_nested_isotopologue_molecules(self, legacy, parent_molecule):
        molecules = []
        for field in ("ice_isos", "ppd_isos", "exgal_isos", "exo_isos"):
            context = field.replace("_isos", "")
            for item in self._as_list(legacy.get(field)):
                if not isinstance(item, dict):
                    continue
                molecules.append(
                    self._convert_nested_isotopologue_molecule(
                        legacy, parent_molecule, context, item
                    )
                )
        return molecules

    def _convert_nested_isotopologue_molecule(
        self, legacy, parent_molecule, context, item
    ):
        original_formula = item.get("formula")
        formula = self._normalize_isotopologue_formula(
            original_formula,
            {
                "__legacy_var": legacy.get("__legacy_var"),
                "__line": item.get("__line"),
                "formula": original_formula,
                "name": legacy.get("name"),
            },
            parent_molecule["label"],
        )
        lab_refs = self._nested_isotopologue_lab_refs(legacy, context, item)

        return {
            "name": parent_molecule["name"],
            "formula": formula,
            "table_formula": self._table_formula(item, original_formula, formula),
            "label": self._isotopologue_label(original_formula),
            "note": self._none_if_blank(item.get("notes")),
            "iupac_name": None,
            "selfies": None,
            "synonyms": [],
            "smiles": self._none_if_blank(item.get("smiles")),
            "canonical_smiles": None,
            "inchi": None,
            "inchikey": None,
            "radical_override": None,
            "fullerene": parent_molecule["fullerene"],
            "pah": parent_molecule["pah"],
            "n_rings": parent_molecule["n_rings"],
            "cyclic": parent_molecule["cyclic"],
            "rotcon": None,
            "dipole": None,
            "refs": {
                "lab": lab_refs,
                "computation": [],
            },
            "isotopologue_of": parent_molecule["label"],
            "latex_header": None,
            "latex_notes": None,
        }

    def _nested_isotopologue_lab_refs(self, legacy, context, item):
        nested_context = dict(item)
        nested_context["__legacy_var"] = legacy.get("__legacy_var")
        nested_context["name"] = legacy.get("name")
        refs = []
        refs.extend(
            self._refs_from_fields(nested_context, "l_ref_bib_ids", "l_refs", "lab")
        )
        refs.extend(
            self._refs_from_fields(
                nested_context,
                f"{context}_l_bib_ids",
                f"{context}_l_refs",
                "lab",
            )
        )
        return self._unique_preserve_order(refs)

    def _isotopologue_label(self, formula):
        text = str(formula or "")
        text = text.replace("^{", "").replace("}", "")
        text = text.replace("[", "").replace("]", "")
        text = re.sub(r"[^A-Za-z0-9+_-]", "", text)
        return f"mol:{text}"

    def _molecule_label(self, legacy):
        raw_label = (
            legacy.get("label")
            or legacy.get("astromol_name")
            or legacy.get("formula")
            or legacy.get("__legacy_var")
        )
        if raw_label is None:
            self._issue(
                "missing_label",
                "error",
                legacy,
                "Legacy molecule has no label, astromol_name, formula, or variable name.",
            )
            raw_label = f"line-{legacy.get('__line')}"
        return raw_label if str(raw_label).startswith("mol:") else f"mol:{raw_label}"

    def _normalize_formula(self, formula, legacy, label):
        normalized = FORMULA_OVERRIDES.get(formula, formula)
        return self._validate_formula(normalized, formula, legacy, label)

    def _normalize_isotopologue_formula(self, formula, legacy, label):
        normalized = ISOTOPOLOGUE_FORMULA_OVERRIDES.get(
            formula, FORMULA_OVERRIDES.get(formula, formula)
        )
        return self._validate_formula(normalized, formula, legacy, label)

    def _validate_formula(self, formula, original_formula, legacy, label):
        try:
            Formula(
                formula,
                parse_groups=False,
                parse_oligos=False,
                parse_fractions=False,
                parse_arithmetic=True,
                allow_empty=False,
            )
        except Exception as exc:
            self._issue(
                "invalid_formula",
                "error",
                legacy,
                f"Formula '{formula}' could not be parsed by molmass: {exc}",
                {
                    "label": label,
                    "formula": formula,
                    "original_formula": original_formula,
                },
            )
        return formula

    def _table_formula(self, legacy, original_formula, formula):
        table_formula = self._none_if_blank(legacy.get("table_formula"))
        if table_formula is not None:
            return table_formula
        if (
            original_formula in FORMULA_OVERRIDES
            or original_formula in ISOTOPOLOGUE_FORMULA_OVERRIDES
        ):
            return original_formula
        return formula

    def _molecule_name(self, legacy, formula, label):
        name = self._none_if_blank(legacy.get("name"))
        if name is not None:
            return name

        name = NAME_OVERRIDES.get(legacy.get("__legacy_var"))
        if name is not None:
            return name

        self._issue(
            "missing_name_filled_from_formula",
            "warning",
            legacy,
            f"Molecule name is blank; using formula '{formula}' in preview.",
            {"label": label, "formula": formula},
        )
        return formula

    def _rotcon(self, legacy, refs):
        rotcon = {
            "A": legacy.get("Acon"),
            "B": legacy.get("Bcon"),
            "C": legacy.get("Ccon"),
            "refs": refs,
            "note": None,
        }
        return rotcon if any(rotcon[key] is not None for key in ("A", "B", "C")) else None

    def _dipole(self, legacy, refs):
        dipole = {
            "a": legacy.get("mua"),
            "b": legacy.get("mub"),
            "c": legacy.get("muc"),
            "refs": refs,
            "note": None,
        }
        return dipole if any(dipole[key] is not None for key in ("a", "b", "c")) else None

    def _convert_main_detection(self, legacy, label):
        sources = self._source_refs(legacy.get("sources"), legacy, "sources")
        telescopes = self._telescope_refs(legacy.get("telescopes"), legacy, "telescopes")
        wavelengths = self._wavelengths(legacy.get("wavelengths"), legacy)
        year = legacy.get("year")

        if year is None:
            self._issue(
                "main_detection_omitted",
                "error",
                legacy,
                "Main ISM/CSM detection has no year.",
                {"label": label},
            )
            return None
        if not sources or not telescopes or not wavelengths:
            self._issue(
                "main_detection_omitted",
                "error",
                legacy,
                "Main ISM/CSM detection is missing required source, telescope, or wavelength data.",
                {"label": label},
            )
            return None

        observation_refs = self._refs_from_fields(
            legacy, "d_ref_bib_ids", "d_refs", "observation"
        )

        return {
            "note": None,
            "molecule": label,
            "sources": sources,
            "telescopes": telescopes,
            "wavelengths": wavelengths,
            "year": year,
            "type": "ISM/CSM",
            "first": True,
            "refs": {
                "observation": observation_refs,
            },
            "latex_text": None,
        }

    def _convert_nested_isotopologue_detections(self, legacy):
        detections = []
        for field in NESTED_ISOTOPOLOGUE_FIELDS:
            context = field.replace("_isos", "")
            for item in self._as_list(legacy.get(field)):
                if not isinstance(item, dict):
                    continue

                detection = self._convert_nested_isotopologue_detection(
                    legacy, field, context, item
                )
                if detection is not None:
                    detections.append(detection)

        return detections

    def _convert_nested_isotopologue_detection(self, legacy, field, context, item):
        nested_context = self._nested_isotopologue_context(legacy, item)
        label = self._isotopologue_label(item.get("formula"))

        source_value = self._context_field_value(legacy, item, context, "sources")
        telescope_value = self._context_field_value(legacy, item, context, "telescopes")
        wavelength_value = self._context_field_value(legacy, item, context, "wavelengths")

        sources = self._source_refs(source_value, nested_context, f"{context}_sources")
        telescopes = self._telescope_refs(
            telescope_value, nested_context, f"{context}_telescopes"
        )
        wavelengths = self._wavelengths(wavelength_value, nested_context)
        observation_refs = self._nested_isotopologue_observation_refs(
            nested_context, context
        )
        year = self._nested_detection_year(legacy, item, context, observation_refs)

        missing_metadata = []
        if not sources:
            missing_metadata.append("sources")
        if not telescopes:
            missing_metadata.append("telescopes")
        if not wavelengths:
            missing_metadata.append("wavelengths")

        missing_required = []
        if year is None:
            missing_required.append("year")
        if not observation_refs:
            missing_required.append("observation_refs")

        detection_metadata = {
            "source_field": source_value,
            "sources": sources,
            "telescope_field": telescope_value,
            "telescopes": telescopes,
            "wavelength_field": wavelength_value,
            "wavelengths": wavelengths,
            "year": year,
            "observation_refs": observation_refs,
        }

        if missing_metadata:
            self._report_nested_isotopologue_detection_metadata_missing(
                legacy,
                item,
                field,
                context,
                label,
                missing_metadata,
                detection_metadata,
            )

        if missing_required:
            self._report_incomplete_nested_isotopologue_detection(
                legacy,
                item,
                field,
                context,
                label,
                missing_required,
                detection_metadata,
            )
            return None

        return {
            "note": self._none_if_blank(item.get("notes") or item.get("note")),
            "molecule": label,
            "sources": sources,
            "telescopes": telescopes,
            "wavelengths": wavelengths,
            "year": year,
            "type": context,
            "first": True,
            "refs": {
                "observation": observation_refs,
            },
            "latex_text": None,
        }

    def _nested_isotopologue_context(self, legacy, item):
        context = dict(item)
        context["__legacy_var"] = legacy.get("__legacy_var")
        context["__line"] = item.get("__line")
        context["formula"] = item.get("formula")
        context["name"] = legacy.get("name")
        return context

    def _context_field_value(self, legacy, item, context, suffix):
        for key in (f"{context}_{suffix}", suffix):
            if key in item:
                return item.get(key)

        key = f"{context}_{suffix}"
        if key in legacy:
            return legacy.get(key)

        return None

    def _nested_isotopologue_observation_refs(self, nested_context, context):
        refs = []
        for id_field, text_field in (
            ("d_ref_bib_ids", "d_refs"),
            (f"{context}_d_bib_ids", f"{context}_d_refs"),
        ):
            refs.extend(
                self._refs_from_fields(nested_context, id_field, text_field, "observation")
            )
        return self._unique_preserve_order(refs)

    def _nested_detection_year(self, legacy, item, context, observation_refs):
        for value in (
            item.get("year"),
            item.get(f"{context}_year"),
            legacy.get(f"{context}_year"),
        ):
            year = self._year_value(value)
            if year is not None:
                return year

        years = [
            self.reference_years[ref]
            for ref in observation_refs
            if ref in self.reference_years
        ]
        return min(years) if years else None

    def _year_value(self, value):
        if value is None:
            return None
        match = re.search(r"\d{4}", str(value))
        return int(match.group(0)) if match else None

    def _report_incomplete_nested_isotopologue_detection(
        self, legacy, item, field, context, label, missing, found
    ):
        self._report_nested_isotopologue_detection_issue(
            "nested_isotopologue_detection_incomplete",
            "Nested isotopologue detection lacks required detection metadata.",
            legacy,
            item,
            field,
            context,
            label,
            missing,
            found,
        )

    def _report_nested_isotopologue_detection_metadata_missing(
        self, legacy, item, field, context, label, missing, found
    ):
        self._report_nested_isotopologue_detection_issue(
            "nested_isotopologue_detection_metadata_missing",
            "Nested isotopologue detection emitted with placeholder metadata.",
            legacy,
            item,
            field,
            context,
            label,
            missing,
            found,
        )

    def _report_nested_isotopologue_detection_issue(
        self, kind, message, legacy, item, field, context, label, missing, found
    ):
        formula = item.get("formula")
        preview_formula = ISOTOPOLOGUE_FORMULA_OVERRIDES.get(
            formula, FORMULA_OVERRIDES.get(formula, formula)
        )
        self._issue(
            kind,
            "info",
            self._nested_isotopologue_context(legacy, item),
            message,
            {
                "parent_line": legacy.get("__line"),
                "parent_formula": legacy.get("formula"),
                "parent_name": legacy.get("name"),
                "parent_label": self._molecule_label(legacy),
                "field": field,
                "context": context,
                "molecule_label": label,
                "formula": formula,
                "preview_formula": preview_formula,
                "table_formula": self._table_formula(item, formula, preview_formula),
                "missing_fields": missing,
                **found,
            },
        )

    def _source_refs(self, values, legacy, field):
        refs = []
        for value in self._as_field_list(values):
            mapped = SOURCE_ALIASES.get(value, value)
            if mapped not in self.source_nicks:
                self._issue(
                    "unknown_source",
                    "error",
                    legacy,
                    f"Unknown source '{value}' in field '{field}'.",
                    {"field": field, "source": value, "mapped_source": mapped},
                )
                continue
            refs.append(mapped)
        return refs

    def _telescope_refs(self, values, legacy, field):
        refs = []
        for value in self._as_field_list(values):
            if value not in self.telescope_nicks:
                self._issue(
                    "unknown_telescope",
                    "error",
                    legacy,
                    f"Unknown telescope '{value}' in field '{field}'.",
                    {"field": field, "telescope": value},
                )
                continue
            refs.append(value)
        return refs

    def _wavelengths(self, values, legacy):
        wavelengths = []
        for value in self._as_field_list(values):
            if value not in ALLOWED_WAVELENGTHS:
                self._issue(
                    "unknown_wavelength",
                    "error",
                    legacy,
                    f"Unknown wavelength '{value}'.",
                    {"wavelength": value},
                )
                continue
            wavelengths.append(value)
        return wavelengths

    def _refs_from_fields(self, legacy, id_field, text_field, role):
        refs = []
        raw_ids = self._as_list(legacy.get(id_field))
        text_chunks = self._split_free_text_refs(legacy.get(text_field))

        if raw_ids:
            for index, ref_id in enumerate(raw_ids):
                ref = self._resolve_reference_key(
                    ref_id, legacy, role, emit_issue=False
                )
                if ref is None and index < len(text_chunks):
                    ref = self._resolve_free_text_reference(
                        text_chunks[index],
                        legacy,
                        role,
                        source_ref=ref_id,
                        emit_issue=False,
                    )
                if ref is None:
                    ref = self._resolve_reference_key(
                        ref_id, legacy, role, emit_issue=True
                    )
                if ref is not None:
                    refs.extend(self._as_list(ref))
            return self._unique_preserve_order(refs)

        for chunk in text_chunks:
            ref = self._resolve_free_text_reference(chunk, legacy, role)
            if ref is not None:
                refs.extend(self._as_list(ref))

        return self._unique_preserve_order(refs)

    def _resolve_reference_key(self, ref_id, legacy, role, emit_issue=True):
        alias = self.reference_key_aliases.get(self._reference_alias_key(ref_id))
        if alias is not None:
            mapped_refs = self._known_alias_refs(alias, ref_id, legacy, role)
            if mapped_refs:
                return mapped_refs

        if ref_id in self.reference_ids:
            return ref_id

        candidates = self._candidate_refs_for_key(ref_id)
        if len(candidates) == 1:
            mapped_ref = candidates[0]
            self._record_reference_mapping(
                ref_id,
                mapped_ref,
                "legacy_key_metadata",
                legacy,
                role,
            )
            return mapped_ref

        if emit_issue:
            if candidates:
                self._issue(
                    "ambiguous_reference",
                    "warning",
                    legacy,
                    f"Reference key '{ref_id}' matches multiple references.bib keys.",
                    {
                        "role": role,
                        "bibcode": ref_id,
                        "candidates": candidates,
                    },
                )
            else:
                self._issue(
                    "unresolved_reference",
                    "warning",
                    legacy,
                    f"Reference key '{ref_id}' does not exist in references.bib.",
                    {"role": role, "bibcode": ref_id},
                )
        return None

    def _candidate_refs_for_key(self, ref_id):
        candidates = []

        base_candidates = self.reference_indexes["zotero_base"].get(ref_id, [])
        candidates.extend(base_candidates)

        ads_parts = self._parse_ads_bibcode(ref_id)
        if ads_parts is not None:
            year, initial, volume, pages = ads_parts
            for page in pages:
                volume_candidates = []
                if volume is not None:
                    volume_candidates = self.reference_indexes["ads_volume"].get(
                        (year, initial, volume, page), []
                    )
                if volume_candidates:
                    candidates.extend(volume_candidates)
                else:
                    candidates.extend(
                        self.reference_indexes["ads"].get((year, initial, page), [])
                    )

        return self._unique_preserve_order(candidates)

    def _parse_ads_bibcode(self, ref_id):
        text = str(ref_id)
        if len(text) != 19 or not text[:4].isdigit():
            return None

        year = text[:4]
        volume = self._normalize_volume(text[9:13])
        qualifier = text[13].replace(".", "")
        page = text[14:18].replace(".", "")
        initial = text[18].upper()
        if not page:
            return None

        pages = {page.upper()}
        if qualifier:
            pages.add(f"{qualifier}{page}".upper())
            pages.add(f"{page}{qualifier}".upper())

        return year, initial, volume, pages

    def _resolve_free_text_reference(
        self, chunk, legacy, role, source_ref=None, emit_issue=True
    ):
        alias = self.free_text_reference_aliases.get(self._reference_alias_key(chunk))
        if alias is not None:
            mapped_refs = self._known_alias_refs(
                alias,
                source_ref or chunk,
                legacy,
                role,
                strategy="free_text_alias",
            )
            if mapped_refs:
                return mapped_refs

        parsed = self._parse_free_text_reference(chunk)
        if parsed is None:
            if emit_issue:
                self._issue(
                    "free_text_reference_unparsed",
                    "warning",
                    legacy,
                    f"Could not parse free-text reference: {chunk}",
                    {
                        "role": role,
                        "source_ref": source_ref,
                        "text": chunk,
                    },
                )
            return None

        year, surname, volume, page = parsed
        candidates = []
        if volume is not None:
            candidates = self.reference_indexes["text_volume"].get(parsed, [])
        if not candidates:
            candidates = self.reference_indexes["text"].get(
                (year, surname, page), []
            )
        candidates = self._unique_preserve_order(candidates)

        if len(candidates) == 1:
            mapped_ref = candidates[0]
            self._record_reference_mapping(
                source_ref or chunk,
                mapped_ref,
                "free_text_metadata",
                legacy,
                role,
            )
            return mapped_ref

        if emit_issue:
            kind = (
                "ambiguous_free_text_reference"
                if candidates
                else "unresolved_free_text_reference"
            )
            message = (
                f"Free-text reference matches multiple references.bib keys: {chunk}"
                if candidates
                else f"Free-text reference does not match references.bib: {chunk}"
            )
            self._issue(
                kind,
                "warning",
                legacy,
                message,
                {
                    "role": role,
                    "source_ref": source_ref,
                    "text": chunk,
                    "parsed": {
                        "year": year,
                        "surname": surname,
                        "volume": volume,
                        "page": page,
                    },
                    "candidates": candidates,
                },
            )
        return None

    def _known_alias_refs(
        self, alias, original, legacy, role, strategy="legacy_key_alias"
    ):
        mapped_refs = []
        missing_refs = []
        for mapped_ref in self._as_list(alias):
            if mapped_ref in self.reference_ids:
                mapped_refs.append(mapped_ref)
                self._record_reference_mapping(
                    original,
                    mapped_ref,
                    strategy,
                    legacy,
                    role,
                )
            else:
                missing_refs.append(mapped_ref)

        if missing_refs:
            self._issue(
                "missing_alias_target",
                "warning",
                legacy,
                f"Reference alias for '{original}' points to missing references.bib keys.",
                {
                    "role": role,
                    "original": original,
                    "mapped": missing_refs,
                },
            )

        return mapped_refs

    def _parse_free_text_reference(self, chunk):
        text = str(chunk).strip()
        year_match = re.search(r"\b(18|19|20)\d{2}\b", text)
        if not year_match:
            return None

        author_text = text[: year_match.start()].strip(" ,")
        year = year_match.group(0)
        volume, page = self._volume_page_from_reference_text(
            text[year_match.end() :]
        )
        if not author_text or not page:
            return None

        author_text = re.split(
            r"\s+(?:et\s+al\.?|and|&)(?:\s+|$)", author_text, maxsplit=1
        )[0]
        surname = self._normalize_lookup_text(author_text)
        if not surname:
            return None

        return year, surname, volume, page

    def _volume_page_from_reference_text(self, text):
        text = text.strip()
        range_match = re.search(
            r"(?:,\s*|pp\.?\s+)([A-Za-z]?\d+[A-Za-z]?)(?:\s*[-–]\s*[A-Za-z]?\d+[A-Za-z]?)?\s*$",
            text,
            re.IGNORECASE,
        )
        if range_match:
            tokens = re.findall(r"\b\d+[A-Za-z]?\b", text[: range_match.start()])
            volume = self._normalize_volume(tokens[-1]) if tokens else None
            page = self._normalize_reference_page(range_match.group(1))
            return volume, page

        tokens = re.findall(r"\b[A-Za-z]?\d+[A-Za-z]?\b", text)
        if not tokens:
            return None, None
        volume = self._normalize_volume(tokens[-2]) if len(tokens) > 1 else None
        page = self._normalize_reference_page(tokens[-1])
        return volume, page

    def _normalize_reference_page(self, page):
        return re.sub(r"[^A-Za-z0-9]", "", str(page)).upper()

    def _split_free_text_refs(self, text):
        text = self._none_if_blank(text)
        if text is None:
            return []
        return [chunk.strip() for chunk in str(text).split(";") if chunk.strip()]

    def _record_reference_mapping(self, original, mapped, strategy, legacy, role):
        self.reference_mappings.append(
            {
                "legacy_var": legacy.get("__legacy_var"),
                "formula": legacy.get("formula"),
                "line": legacy.get("__line"),
                "role": role,
                "original": original,
                "mapped": mapped,
                "strategy": strategy,
            }
        )

    def _unique_preserve_order(self, values):
        seen = set()
        unique = []
        for value in values:
            if value not in seen:
                seen.add(value)
                unique.append(value)
        return unique

    def _report_omitted_extra_contexts(self, legacy):
        for context in EXTRA_CONTEXTS:
            flag = legacy.get(context)
            if flag not in (True, "Tentative"):
                continue

            self._issue(
                "extra_context_detection_omitted",
                "info",
                legacy,
                f"Legacy '{context}' detection omitted from preview.",
                {
                    "context": context,
                    "flag": flag,
                    "source_field": legacy.get(f"{context}_sources"),
                    "telescope_field": legacy.get(f"{context}_telescopes"),
                    "wavelength_field": legacy.get(f"{context}_wavelengths"),
                    "refs": legacy.get(f"{context}_d_bib_ids"),
                    "free_text_refs": legacy.get(f"{context}_d_refs"),
                },
            )

            self._refs_from_fields(
                legacy,
                f"{context}_d_bib_ids",
                f"{context}_d_refs",
                f"{context}_observation",
            )
            self._refs_from_fields(
                legacy,
                f"{context}_l_bib_ids",
                f"{context}_l_refs",
                f"{context}_lab",
            )

    def _validate_preview(self, molecules, detections):
        labels = [molecule["label"] for molecule in molecules]
        for label, count in Counter(labels).items():
            if count > 1:
                self._issue(
                    "duplicate_preview_label",
                    "error",
                    {"__legacy_var": label, "__line": None},
                    f"Preview label '{label}' appears {count} times.",
                )

        molecule_labels = set(labels)
        for detection in detections:
            if detection["molecule"] not in molecule_labels:
                self._issue(
                    "preview_detection_unknown_molecule",
                    "error",
                    {"__legacy_var": detection["molecule"], "__line": None},
                    f"Preview detection references unknown molecule '{detection['molecule']}'.",
                )

    def _issue(self, kind, severity, legacy, message, details=None):
        self.issues.append(
            {
                "kind": kind,
                "severity": severity,
                "legacy_var": legacy.get("__legacy_var"),
                "name": legacy.get("name"),
                "formula": legacy.get("formula"),
                "line": legacy.get("__line"),
                "message": message,
                "details": details or {},
            }
        )

    def _issue_report(self, molecules, detections):
        by_kind = Counter(issue["kind"] for issue in self.issues)
        by_severity = Counter(issue["severity"] for issue in self.issues)
        mappings_by_strategy = Counter(
            mapping["strategy"] for mapping in self.reference_mappings
        )
        return {
            "summary": {
                "converted_molecules": len(molecules),
                "converted_detections": len(detections),
                "issues": len(self.issues),
                "issues_by_kind": dict(sorted(by_kind.items())),
                "issues_by_severity": dict(sorted(by_severity.items())),
                "reference_mappings": len(self.reference_mappings),
                "reference_mappings_by_strategy": dict(
                    sorted(mappings_by_strategy.items())
                ),
            },
            "reference_mappings": self.reference_mappings,
            "issues": self.issues,
        }

    def _as_list(self, value):
        if value is None:
            return []
        if isinstance(value, list):
            return value
        return [value]

    def _as_field_list(self, value):
        if isinstance(value, str):
            return [part.strip() for part in value.split(",") if part.strip()]
        return self._as_list(value)

    def _none_if_blank(self, value):
        if value == "":
            return None
        return value


def write_json(path, data):
    path.write_text(json.dumps(data, indent=2, ensure_ascii=True) + "\n")


def triage_value(value):
    if value is None:
        return None
    if isinstance(value, list):
        return ", ".join(f"`{item}`" for item in value) if value else "[]"
    if isinstance(value, dict):
        return f"`{json.dumps(value, ensure_ascii=False, sort_keys=True)}`"
    return f"`{value}`"


def write_triage(path, report):
    issues = report["issues"]
    summary = report["summary"]
    mappings = report["reference_mappings"]

    lines = [
        "# Legacy Conversion Triage",
        "",
        "Generated from `scripts/convert_legacy_molecules.py`.",
        "",
        "## Summary",
        "",
        f"- Converted molecules: {summary['converted_molecules']}",
        f"- Converted detections: {summary['converted_detections']}",
        f"- Reference mappings applied: {summary['reference_mappings']}",
        f"- Remaining issues: {summary['issues']}",
        "",
        "## Remaining Issues by Kind",
        "",
    ]

    for kind, count in summary["issues_by_kind"].items():
        lines.append(f"- `{kind}`: {count}")

    lines.extend(["", "## Reference Alias Mappings", ""])
    alias_mappings = [
        mapping
        for mapping in mappings
        if mapping["strategy"] in {"legacy_key_alias", "free_text_alias"}
    ]
    if alias_mappings:
        for mapping in alias_mappings:
            lines.append(
                "- "
                f"line {mapping['line']} `{mapping['legacy_var']}` "
                f"[{mapping['role']}] {mapping['original']} -> {mapping['mapped']}"
            )
    else:
        lines.append("- None")

    lines.extend(["", "## Ambiguous References", ""])
    ambiguous = [
        issue
        for issue in issues
        if issue["kind"] in {"ambiguous_reference", "ambiguous_free_text_reference"}
    ]
    if ambiguous:
        for issue in ambiguous:
            details = issue["details"]
            ref = details.get("bibcode") or details.get("text")
            candidates = ", ".join(details.get("candidates", []))
            lines.append(
                f"- line {issue['line']} `{issue['legacy_var']}`: {ref} -> {candidates}"
            )
    else:
        lines.append("- None")

    lines.extend(["", "## Unresolved Reference Keys", ""])
    unresolved_refs = defaultdict(list)
    for issue in issues:
        if issue["kind"] == "unresolved_reference":
            unresolved_refs[issue["details"]["bibcode"]].append(issue)

    if unresolved_refs:
        for ref, grouped_issues in sorted(unresolved_refs.items()):
            lines.append(f"- {len(grouped_issues)}x `{ref}`")
            for issue in grouped_issues:
                details = issue["details"]
                name = f" name `{issue['name']}`" if issue.get("name") else ""
                lines.append(
                    "  - "
                    f"line {issue['line']} `{issue['legacy_var']}` "
                    f"formula `{issue['formula']}`{name} [{details.get('role')}]"
                )
    else:
        lines.append("- None")

    lines.extend(["", "## Unresolved Free-Text References", ""])
    unresolved_text = defaultdict(list)
    for issue in issues:
        if issue["kind"] == "unresolved_free_text_reference":
            unresolved_text[issue["details"]["text"]].append(issue)

    if unresolved_text:
        for text, grouped_issues in sorted(unresolved_text.items()):
            lines.append(f"- {len(grouped_issues)}x {text}")
            for issue in grouped_issues:
                details = issue["details"]
                source_ref = details.get("source_ref")
                source_text = f", source `{source_ref}`" if source_ref else ""
                name = f" name `{issue['name']}`" if issue.get("name") else ""
                lines.append(
                    "  - "
                    f"line {issue['line']} `{issue['legacy_var']}` "
                    f"formula `{issue['formula']}`{name} [{details.get('role')}]"
                    f"{source_text}"
                )
    else:
        lines.append("- None")

    lines.extend(["", "## Unknown Data References", ""])
    unknown_issues = [
        issue
        for issue in issues
        if issue["kind"] in {"unknown_source", "unknown_telescope", "unknown_wavelength"}
    ]
    if unknown_issues:
        for issue in unknown_issues:
            details = issue["details"]
            value = (
                details.get("source")
                or details.get("telescope")
                or details.get("wavelength")
            )
            mapped = details.get("mapped_source")
            mapped_text = f"; mapped `{mapped}`" if mapped and mapped != value else ""
            lines.append(
                "- "
                f"line {issue['line']} `{issue['legacy_var']}` "
                f"formula `{issue['formula']}`; field `{details.get('field')}`; "
                f"value `{value}`{mapped_text}"
            )
    else:
        lines.append("- None")

    lines.extend(["", "## Nested Isotopologue Detection Gaps", ""])
    nested_issues = [
        issue
        for issue in issues
        if issue["kind"]
        in {
            "nested_isotopologue_detection_incomplete",
            "nested_isotopologue_detection_metadata_missing",
        }
    ]
    if nested_issues:
        for issue in nested_issues:
            details = issue["details"]
            parent_name = (
                f" name `{details['parent_name']}`"
                if details.get("parent_name")
                else ""
            )
            missing = triage_value(details.get("missing_fields"))
            lines.append(
                "- "
                f"parent line {details['parent_line']} `{issue['legacy_var']}` "
                f"formula `{details['parent_formula']}`{parent_name}; "
                f"item line {issue['line']}; label `{details['molecule_label']}`; "
                f"formula `{details['formula']}`; preview `{details['preview_formula']}`; "
                f"table `{details['table_formula']}`; context `{details['context']}`; "
                f"missing {missing}"
            )

            fields = [
                f"`sources`={triage_value(details.get('source_field'))} -> "
                f"{triage_value(details.get('sources'))}",
                f"`telescopes`={triage_value(details.get('telescope_field'))} -> "
                f"{triage_value(details.get('telescopes'))}",
                f"`wavelengths`={triage_value(details.get('wavelength_field'))} -> "
                f"{triage_value(details.get('wavelengths'))}",
                f"`year`=`{details.get('year')}`",
            ]
            lines.append(f"  - fields: {'; '.join(fields)}")
            lines.append(
                "  - refs: "
                f"resolved {triage_value(details.get('observation_refs'))}"
            )
    else:
        lines.append("- None")

    lines.extend(["", "## Deferred Categories", ""])
    for kind in (
        "extra_context_detection_omitted",
        "nested_isotopologue_detection_incomplete",
        "nested_isotopologue_detection_metadata_missing",
        "missing_name_filled_from_formula",
        "formula_normalized",
    ):
        count = sum(1 for issue in issues if issue["kind"] == kind)
        if count:
            lines.append(f"- `{kind}`: {count}")

    path.write_text("\n".join(lines) + "\n")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--only",
        nargs="*",
        help="Optional legacy variable, label, astromol_name, formula, or mol: label filter.",
    )
    parser.add_argument("--molecules-out", type=Path, default=DEFAULT_MOLECULES_OUT)
    parser.add_argument("--detections-out", type=Path, default=DEFAULT_DETECTIONS_OUT)
    parser.add_argument("--issues-out", type=Path, default=DEFAULT_ISSUES_OUT)
    parser.add_argument("--triage-out", type=Path, default=DEFAULT_TRIAGE_OUT)
    return parser.parse_args()


def main():
    args = parse_args()
    converter = LegacyConverter(only=set(args.only) if args.only else None)
    molecules, detections, issues = converter.convert()

    write_json(args.molecules_out, molecules)
    write_json(args.detections_out, detections)
    write_json(args.issues_out, issues)
    write_triage(args.triage_out, issues)

    print(f"Wrote {len(molecules)} molecules to {args.molecules_out}")
    print(f"Wrote {len(detections)} detections to {args.detections_out}")
    print(f"Wrote {issues['summary']['issues']} issues to {args.issues_out}")
    print(f"Wrote triage report to {args.triage_out}")
    print(json.dumps(issues["summary"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
