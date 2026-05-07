import json
import re
from pathlib import Path
from .models import (
    Ref, Telescope, Source,
    RotationalConstants, DipoleMoment, Molecule, Detection,
    RecordHistory,
    DETECTION_RELATION_FIELDS,
)

# Path to the data directory, relative to this file
DATA_DIR = Path(__file__).parent / "data"

MONTHS = {
    "jan": 1,
    "january": 1,
    "feb": 2,
    "february": 2,
    "mar": 3,
    "march": 3,
    "apr": 4,
    "april": 4,
    "may": 5,
    "jun": 6,
    "june": 6,
    "jul": 7,
    "july": 7,
    "aug": 8,
    "august": 8,
    "sep": 9,
    "sept": 9,
    "september": 9,
    "oct": 10,
    "october": 10,
    "nov": 11,
    "november": 11,
    "dec": 12,
    "december": 12,
}

MOLMASS_DERIVED_FIELDS = {
    "atom_counts",
    "atoms",
    "isotope_counts",
    "mass",
    "average_mass",
    "nominal_mass",
    "natoms",
    "charge",
    "cation",
    "anion",
    "neutral",
    "nelectrons",
    "odd_electron",
}

FORMULA_HEURISTIC_FIELDS = {
    "du",
    "maxdu",
}

SPECTROSCOPY_DERIVED_FIELDS = {
    "is_linear",
    "kappa",
}

MOLECULE_NON_INIT_FIELDS = (
    MOLMASS_DERIVED_FIELDS
    | FORMULA_HEURISTIC_FIELDS
    | SPECTROSCOPY_DERIVED_FIELDS
)

ROTCON_ALIASES = {
    "Acon": "A",
    "Bcon": "B",
    "Ccon": "C",
}

DIPOLE_ALIASES = {
    "mua": "a",
    "mub": "b",
    "muc": "c",
}


class Database:
    """
    Central container for all astromol data.
    Loads JSON files, resolves cross-references, and exposes simple lookups.
    """
    def __init__(self):
        # Storage — these get populated by _load()
        self.refs = {}           # BibTeX citation key -> Ref
        self.telescopes = {}     # nick -> Telescope
        self.sources = {}        # nick -> Source
        self.molecules = {}      # molecule label -> Molecule
        self.molecules_by_formula = {}  # formula -> list[Molecule]
        self.molecules_by_name = {}     # name -> list[Molecule]
        self.detections = []     # list of Detection objects
        self.detections_by_id = {}  # stable detection ID -> Detection

        self._load()

    # ================================================================
    # Loading — called once during __init__
    # ================================================================

    def _load(self):
        """Load all data files in dependency order."""
        self._load_refs()          # no dependencies
        self._load_telescopes()    # no dependencies
        self._load_sources()       # no dependencies
        self._load_molecules()     # depends on refs
        self._load_detections()    # depends on refs, telescopes, sources, molecules

    def _load_refs(self):
        """Load references.bib into self.refs."""
        try:
            import bibtexparser
        except ImportError as exc:
            raise ImportError(
                "Loading references.bib requires the 'bibtexparser' package."
            ) from exc

        with open(DATA_DIR / "references.bib") as f:
            bib_database = bibtexparser.load(f)

        seen = set()
        for entry in bib_database.entries:
            bibcode = entry.get("ID")
            if bibcode in seen:
                raise ValueError(f"Duplicate reference key in references.bib: {bibcode}")
            seen.add(bibcode)

            ref = self._ref_from_bibtex_entry(entry)
            self.refs[ref.bibcode] = ref

    def _ref_from_bibtex_entry(self, entry):
        """Convert a parsed BibTeX entry to a Ref object."""
        bibcode = entry["ID"]
        year = self._parse_int(entry.get("year"))
        if year is None:
            raise ValueError(f"Reference '{bibcode}' is missing a valid year")

        return Ref(
            bibcode=bibcode,
            author=self._format_bibtex_authors(entry),
            journal=self._bibtex_publication(entry),
            year=year,
            month=self._parse_month(entry.get("month")),
            volume=entry.get("volume"),
            page=entry.get("pages") or entry.get("page"),
            title=entry.get("title"),
            doi=entry.get("doi"),
            note=entry.get("annotation"),
        )

    def _format_bibtex_authors(self, entry):
        """Return a compact author string from BibTeX author/editor fields."""
        authors = entry.get("author") or entry.get("editor")
        if not authors:
            return "Unknown"

        names = [self._format_bibtex_name(name) for name in authors.split(" and ")]
        if len(names) == 1:
            return names[0]
        if len(names) == 2:
            return f"{names[0]} & {names[1]}"
        return f"{names[0]} et al."

    def _format_bibtex_name(self, name):
        """Format one BibTeX name as a display surname."""
        name = name.strip().strip("{}")
        if "," in name:
            return name.split(",", 1)[0].strip("{} ")
        return name.split()[-1].strip("{} ")

    def _bibtex_publication(self, entry):
        """Return the best available publication/container string."""
        return (
            entry.get("journal")
            or entry.get("booktitle")
            or entry.get("publisher")
            or entry.get("ENTRYTYPE")
            or "Unknown"
        )

    def _parse_int(self, value):
        """Parse an integer from a BibTeX scalar."""
        if value is None:
            return None
        match = re.search(r"\d+", str(value))
        return int(match.group(0)) if match else None

    def _parse_month(self, value):
        """Parse BibTeX month values into 1-12 integers."""
        if value is None:
            return None
        if isinstance(value, int):
            return value

        text = str(value).strip().strip("{}").lower()
        if text.isdigit():
            return int(text)
        return MONTHS.get(text[:3], MONTHS.get(text))

    def _load_telescopes(self):
        """Load telescopes.json into self.telescopes."""
        with open(DATA_DIR / "telescopes.json") as f:
            for entry in json.load(f):
                entry = {k: v for k, v in entry.items() if not k.startswith("_")}
                self._normalize_history(entry)
                tel = Telescope(**entry)
                self.telescopes[tel.nick] = tel

    def _load_sources(self):
        """Load sources.json into self.sources."""
        with open(DATA_DIR / "sources.json") as f:
            for entry in json.load(f):
                entry = {k: v for k, v in entry.items() if not k.startswith("_")}
                self._normalize_history(entry)
                src = Source(**entry)
                self.sources[src.nick] = src

    def _load_molecules(self):
        """Load molecules.json into label-keyed storage and secondary indexes."""
        with open(DATA_DIR / "molecules.json") as f:
            data = json.load(f)
            if isinstance(data, dict):
                data = [data]

            for entry in data:
                entry = {k: v for k, v in entry.items() if not k.startswith("_")}
                entry = self._normalize_molecule_entry(entry)

                # Build RotationalConstants from nested dict, if present
                rotcon_data = entry.pop("rotcon", None)
                if rotcon_data is not None:
                    rotcon_data = {
                        k: v for k, v in rotcon_data.items() if not k.startswith("_")
                    }
                    refs = self._resolve_spectroscopy_refs(rotcon_data)
                    rotcon = RotationalConstants(**rotcon_data)
                    if refs is not None:
                        rotcon.refs = refs
                    entry["rotcon"] = rotcon

                # Build DipoleMoment from nested dict, if present
                dipole_data = entry.pop("dipole", None)
                if dipole_data is not None:
                    dipole_data = {
                        k: v for k, v in dipole_data.items() if not k.startswith("_")
                    }
                    refs = self._resolve_spectroscopy_refs(dipole_data)
                    dipole = DipoleMoment(**dipole_data)
                    if refs is not None:
                        dipole.refs = refs
                    entry["dipole"] = dipole

                # Resolve molecule-level refs
                # {"lab": ["bib1", "bib2"]} -> {"lab": [Ref, Ref]}
                refs_data = entry.pop("refs", None)
                if refs_data is not None:
                    entry["refs"] = self._resolve_refs_by_role(
                        refs_data,
                        f"molecule '{entry.get('formula')}'",
                    )

                mol = Molecule(**entry)
                if mol.label in self.molecules:
                    raise ValueError(f"Duplicate molecule label: {mol.label}")

                self.molecules[mol.label] = mol
                self.molecules_by_formula.setdefault(mol.formula, []).append(mol)
                self.molecules_by_name.setdefault(mol.name, []).append(mol)

    def _resolve_spectroscopy_refs(self, data):
        """
        Resolve nested spectroscopy references.

        Supports the current ``refs`` list plus legacy ``ref`` values. A legacy
        ``ref`` may be either a single BibTeX key or a list of BibTeX keys.
        """
        refs_data = data.pop("refs", None)
        ref_data = data.pop("ref", None)

        if refs_data is None:
            refs_data = ref_data
        elif ref_data is not None:
            refs_data = self._as_list(refs_data) + self._as_list(ref_data)

        if refs_data is None:
            return None

        refs = []
        for bibcode in self._as_list(refs_data):
            try:
                refs.append(self.refs[bibcode])
            except KeyError as exc:
                raise KeyError(
                    f"Unknown spectroscopy reference key: {bibcode}"
                ) from exc

        return refs

    def _resolve_refs_by_role(self, refs_data, context):
        """Resolve role-keyed BibTeX reference keys to Ref objects."""
        resolved_refs = {}

        for role, bibcodes in refs_data.items():
            if role.startswith("_"):
                continue

            resolved_refs[role] = []
            for bibcode in self._as_list(bibcodes):
                try:
                    resolved_refs[role].append(self.refs[bibcode])
                except KeyError as exc:
                    raise KeyError(
                        f"Unknown {context} reference key: {bibcode}"
                    ) from exc

        return resolved_refs

    def _as_list(self, value):
        """Normalize scalar or list reference data to a list."""
        if value is None:
            return []
        if isinstance(value, list):
            return value
        return [value]

    def _normalize_molecule_entry(self, entry):
        """
        Normalize molecule JSON into the current Molecule dataclass shape.

        Older data may still contain fields that are now Molecule properties,
        or the former flat spectroscopy keys. Keep curated values where they
        still map to explicit dataclass fields and drop values that are not
        accepted by Molecule.__init__.
        """
        entry = dict(entry)
        self._normalize_history(entry)

        if "radical" in entry:
            if "radical_override" not in entry:
                entry["radical_override"] = entry["radical"]
            entry.pop("radical")

        if "latex_header" in entry and "latex_section_override" not in entry:
            entry["latex_section_override"] = entry["latex_header"]
        entry.pop("latex_header", None)

        if "latex_notes" in entry and "latex_body" not in entry:
            entry["latex_body"] = entry["latex_notes"]
        entry.pop("latex_notes", None)

        for key in MOLECULE_NON_INIT_FIELDS:
            entry.pop(key, None)

        self._move_legacy_nested_fields(entry, "rotcon", ROTCON_ALIASES)
        self._move_legacy_nested_fields(entry, "dipole", DIPOLE_ALIASES)

        return entry

    def _move_legacy_nested_fields(self, entry, target_key, aliases):
        """
        Move legacy flat molecule fields into a nested dataclass payload.

        If the current nested object is already present, it wins. Any legacy
        aliases are removed so Molecule does not receive stale constructor
        arguments.
        """
        nested = entry.get(target_key)
        legacy_values = {}

        for old_key, new_key in aliases.items():
            if old_key in entry:
                legacy_values[new_key] = entry.pop(old_key)

        if (
            nested is None
            and any(value is not None for value in legacy_values.values())
        ):
            entry[target_key] = legacy_values

    def _load_detections(self):
        """Load detections.json into self.detections."""
        with open(DATA_DIR / "detections.json") as f:
            for entry in json.load(f):
                entry = {k: v for k, v in entry.items() if not k.startswith("_")}
                self._normalize_history(entry)

                # Resolve molecule label to Molecule object
                entry["molecule"] = self.molecules[entry["molecule"]]

                # Resolve source nicks to Source objects
                entry["sources"] = [
                    self.sources[nick] for nick in entry["sources"]
                ]

                # Resolve telescope nicks to Telescope objects
                entry["telescopes"] = [
                    self.telescopes[nick] for nick in entry["telescopes"]
                ]

                # Resolve BibTeX keys to Ref objects
                refs_data = entry.pop("refs", None)
                if refs_data is not None:
                    entry["refs"] = self._resolve_refs_by_role(
                        refs_data,
                        f"detection of '{entry.get('molecule')}'",
                    )

                det = Detection(**entry)
                if det.id in self.detections_by_id:
                    raise ValueError(f"Duplicate detection id in detections.json: {det.id}")
                self.detections.append(det)
                self.detections_by_id[det.id] = det

        self._validate_detection_relationships()

    def _validate_detection_relationships(self):
        """Ensure detection relationship IDs point to known detections."""
        for det in self.detections:
            for field_name in DETECTION_RELATION_FIELDS:
                for detection_id in getattr(det, field_name):
                    if detection_id not in self.detections_by_id:
                        raise KeyError(
                            f"Detection '{det.id}' has unknown "
                            f"{field_name} target: {detection_id}"
                        )

    def _normalize_history(self, entry):
        """Convert a nested history payload into a RecordHistory object."""
        history_data = entry.get("history")
        if history_data is not None and not isinstance(history_data, RecordHistory):
            history_data = dict(history_data)
            history_data.pop("last_reviewed", None)
            entry["history"] = RecordHistory(**history_data)

    # ================================================================
    # Simple accessors
    # ================================================================

    def get_ref(self, bibcode):
        """Look up a Ref by BibTeX citation key."""
        return self.refs[bibcode]

    def get_telescope(self, nick):
        """Look up a Telescope by nick."""
        return self.telescopes[nick]

    def get_source(self, nick):
        """Look up a Source by nick."""
        return self.sources[nick]

    def get_molecule(self, label):
        """Look up a Molecule by unique label."""
        return self.molecules[label]

    def get_molecules_by_formula(self, formula):
        """Look up all molecules with a given chemical formula."""
        return list(self.molecules_by_formula.get(formula, []))

    def get_molecules_by_name(self, name):
        """Look up all molecules with a given name."""
        return list(self.molecules_by_name.get(name, []))

    def __repr__(self):
        """Summary shown when printing the database."""
        return (
            "astromol Database: "
            f"{len(self.molecules)} molecules, "
            f"{len(self.detections)} detections, "
            f"{len(self.sources)} sources, "
            f"{len(self.telescopes)} telescopes, "
            f"{len(self.refs)} references"
        )
