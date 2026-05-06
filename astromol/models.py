from dataclasses import dataclass, field
from datetime import date
from functools import cached_property
import re
from molmass import ELECTRON, ELEMENTS, Formula, FormulaError

# from rdkit import Chem # type: ignore (VSCode is a terrible program)
# from rdkit.Chem import Descriptors # type: ignore (VSCode is a terrible program)



# ============================================================
# Allowed values — add new entries here as needed
# ============================================================

SOURCE_TYPES = [
    "SFR",
    "Dark Cloud",
    "Carbon Star",
    "Oxygen Star",
    "Planetary Nebula",
    "Protostellar",
    "Diffuse Cloud",
    "HII Region",
    "PDR",
    "Shock",
    "Supernova Remnant",
    "Sgr A",
    "External Galaxy",
    "YSO",
    "Exoplanet",
    "Other",
]

MOLECULE_REF_ROLES = [
    "lab",            # laboratory spectroscopy papers
    "computation",    # theoretical/computational chemistry papers
]

DETECTION_TYPES = [
    "ISM/CSM",        # interstellar/circumstellar medium
    "isotopologue",   # detection of an isotopologue
    "ice",            # ice/solid-phase detection
    "exgal",          # extragalactic detection
    "ppd",            # protoplanetary disk detection
    "exo",            # exoplanet atmosphere detection
    "comet",          # detection in a comet
    "tentative",      # tentative/unconfirmed detection
]

WAVELENGTHS = [
    "cm",             # centimeter
    "mm",             # millimeter
    "sub-mm",         # sub-millimeter
    "IR",             # infrared
    "UV",             # ultraviolet
    "Vis",            # visible
]

SIMPLE_FORMULA_TOKEN = re.compile(
    r"(?:\[(?P<isotope>\d+)(?P<isotope_symbol>[A-Z][a-z]?)\]|"
    r"(?P<symbol>[A-Z][a-z]?))(?P<count>\d*)"
)

ISOTOPE_MASS_OVERRIDES = {
    # Relative atomic masses for isotope labels used by astromol but not
    # included in molmass' isotope table. 26Al is from AME2020 data as exposed
    # by the periodictable package.
    "26Al": 25.98689188,
}

DETECTION_REF_ROLES = [
    "observation",    # papers reporting the astronomical detection
    "confirmation",   # paper confirming a tentative detection
    "refutation",     # paper challenging a detection
    "correction",     # paper correcting/retracting a detection
]


@dataclass
class Ref:
    """A single bibliographic reference."""

    # Required fields
    bibcode: str          # e.g. "Swings:1937:483" — BibTeX citation key
    author: str           # e.g. "Swings & Rosenfeld"
    journal: str          # e.g. "ApJ"
    year: int             # e.g. 1937 — always known

    # Optional date precision — fill in what you know
    month: int = None     # 1-12, if known
    day: int = None       # 1-31, if known

    # Optional metadata
    volume: str = None    # e.g. "86"
    page: str = None     # e.g. "483"
    title: str = None     # full title of the paper
    doi: str = None       # e.g. "10.1086/143882"
    note: str = None          # private notes, from _note in JSON

    @property
    def sortdate(self):
        """
        Build a date for sorting/ordering purposes.
        Uses the 1st for any unknown month or day.
        Two papers from 1937 with no month both get Jan 1 1937,
        but a paper from March 1937 sorts after January 1937.
        """
        m = self.month if self.month is not None else 1
        d = self.day if self.day is not None else 1
        return date(self.year, m, d)

    def __repr__(self):
        """What you see when you print this object."""
        return f"{self.author} {self.year}, {self.journal}"
    
@dataclass
class Telescope:
    """An observing facility or instrument."""

    # Required fields
    name: str             # e.g. "Green Bank Telescope"
    nick: str             # e.g. "GBT" — unique key, used in detections
    shortname: str        # e.g. "GBT 100-m" — human-readable short name
    type: str             # e.g. "Single Dish", "Interferometer", "Space"
    wavelength: list      # e.g. ["cm", "mm"] — operational wavelength regimes

    # Optional fields
    diameter: float = None        # meters (null for arrays like ALMA)
    latitude: float = None        # decimal degrees
    longitude: float = None       # decimal degrees
    built: int = None             # year constructed/commissioned
    decommissioned: int = None    # year decommissioned, null if still active
    note: str = None              # any free-text notes
    latex_name: str = None        # e.g. "GBT" — for paper generation

    @property
    def active(self):
        """True if the telescope has not been decommissioned."""
        return self.decommissioned is None

    def __repr__(self):
        """What you see when you print this object."""
        return self.nick

@dataclass
class Source:
    """An astronomical object where molecules have been detected."""

    # Required fields
    name: str             # e.g. "Sgr B2" — human-readable name
    nick: str             # e.g. "SgrB2" — unique key, used in detections
    type: str             # must be one of SOURCE_TYPES above

    # Optional fields
    ra: str = None        # right ascension, e.g. "17:47:20.0"
    dec: str = None       # declination, e.g. "-28:23:07.0"
    simbad_url: str = None  # link to SIMBAD entry
    latex_name: str = None  # e.g. "Sgr~B2" — for paper generation
    note: str = None          # private notes, from _note in JSON

    def __post_init__(self):
        """Runs automatically after __init__. Validates the data."""
        if self.type not in SOURCE_TYPES:
            raise ValueError(
                f"Source '{self.nick}': unknown type '{self.type}'. "
                f"Must be one of: {SOURCE_TYPES}"
            )

    def __repr__(self):
        """What you see when you print this object."""
        return self.nick

@dataclass
class RotationalConstants:
    """Rotational constants of a molecule in MHz."""

    A: float = None       # A rotational constant (MHz) — largest
    B: float = None       # B rotational constant (MHz) — middle
    C: float = None       # C rotational constant (MHz) — smallest
    refs: list = field(default_factory=list) # BibTeX keys of source papers
                                             # resolved to Ref objects during loading
    note: str = None

    @property
    def ref(self):
        """Legacy alias for the first reference, if present."""
        return self.refs[0] if self.refs else None

    def __repr__(self):
        """What you see when you print this object."""
        if self.A is None and self.C is None and self.B is not None:
            return f"B={self.B} MHz"
        return f"A={self.A}, B={self.B}, C={self.C} MHz"


@dataclass
class DipoleMoment:
    """Dipole moment components of a molecule in Debye."""

    a: float = None       # a-component (Debye)
    b: float = None       # b-component (Debye)
    c: float = None       # c-component (Debye)
    refs: list = field(default_factory=list) # BibTeX keys of source papers
                                             # resolved to Ref objects during loading
    note: str = None

    @property
    def ref(self):
        """Legacy alias for the first reference, if present."""
        return self.refs[0] if self.refs else None

    @property
    def total(self):
        """Total dipole moment magnitude in Debye."""
        components = [x for x in [self.a, self.b, self.c] if x is not None]
        if not components:
            return None
        return sum(x**2 for x in components) ** 0.5

    def __repr__(self):
        """What you see when you print this object."""
        parts = []
        if self.a is not None:
            parts.append(f"μa={self.a}")
        if self.b is not None:
            parts.append(f"μb={self.b}")
        if self.c is not None:
            parts.append(f"μc={self.c}")
        return ", ".join(parts) + " D" if parts else "unknown"

@dataclass
class Molecule:
    """A chemical species detected (or potentially detectable) in the ISM/CSM.

    The optional ``smiles`` and ``selfies`` fields are retained as identifiers
    or future display metadata only; they are stored but not interpreted here.
    """

    # === Identity ===
    name: str                     # e.g. "methylidyne" — common name
    formula: str                  # e.g. "CH" — chemical formula, not necessarily unique
    

    # === Display / LaTeX ===
    table_formula: str = None     # e.g. "\ce{CH}" — mhchem formatted for LaTeX tables
                                  # defaults to formula if not provided
    label: str = None             # e.g. "mol:CH" — REQUIRED, for \label{}/\ref{} in LaTeX
                                  # must be unique, even for isomers with the same formula

    # === Notes ===
    note: str = None              # private notes, from _note in JSON

    # === Optional identifiers ===
    iupac_name: str = None        # e.g. "methylidyne" — formal IUPAC name
    selfies: str = None           # SELFIES string (alternative to SMILES)
    synonyms: list[str] = field(default_factory=list) #list of synonyms
    smiles: str = None            # SMILES string e.g. "[CH]"
    canonical_smiles: str = None  # Canonical SMILES string
    inchi: str = None             # IUPAC International Chemical Identifier
    inchikey: str = None          # A 27-character, fixed-length hashed version of the full IUPAC International Chemical Identifier (InChI)

    # === Flags that require chemical judgment ===
    # A formula alone cannot determine connectivity, topology, PAH/fullerene
    # classification, or complete electronic structure. The radical flag is
    # therefore implemented as an override below: if set, it wins; if left as
    # None, the odd-electron rule is used as a useful default.
    radical_override: bool = None # True/False to override odd-electron inference
    fullerene: bool = False       # e.g. True for C60
    pah: bool = False             # True for polycyclic aromatic hydrocarbons
    n_rings: int = 0              # Number of rings in the molecule
    cyclic: bool = False          # Is the molecule cyclic or not
    # Manually curated analysis tags, grouped by category.
    tags: dict[str, list[str]] = field(default_factory=dict)

    # === Spectroscopic constants ===
    # Lab/computational data — can't be computed from formula.
    # Each has its own references for where the data came from.
    rotcon: RotationalConstants = None   # rotational constants (MHz) - their own object (above)
    dipole: DipoleMoment = None          # dipole moment components (Debye) - their own object (above)

    # === Molecule-level references ===
    # About the species itself, not any particular detection.
    # Keys are roles from MOLECULE_REF_ROLES, values are lists of BibTeX keys.
    # Resolved to lists of Ref objects during loading.
    refs: dict = None             # e.g. {"lab": ["bib1"], "computation": ["bib2"]}

    # === Isotopologue relationship ===
    isotopologue_of: str = None   # label of parent molecule, e.g. "mol:CH" for 13CH
                                  # null means this IS a parent molecule

    # === LaTeX generation ===
    latex_section_override: str = None  # optional custom section heading
    latex_body: str = None              # curated prose for the paper section

    def __post_init__(self):
        """Validate metadata and fill in simple defaults."""

        # Default table_formula to formula if not provided.
        if self.table_formula is None:
            self.table_formula = self.formula

        # label is required — isomers can share a formula so labels must be unique.
        if self.label is None:
            raise ValueError(
                f"Molecule '{self.formula}': label is required. "
                f"Set it to a unique LaTeX label string, e.g. 'mol:CH'"
            )

        # Default refs to empty dict.
        if self.refs is None:
            self.refs = {}

        # Validate manually curated tag metadata.
        if self.tags is None:
            self.tags = {}
        if not isinstance(self.tags, dict):
            raise ValueError(
                f"Molecule '{self.formula}': tags must be a dict of lists."
            )
        for category, values in self.tags.items():
            if not isinstance(category, str):
                raise ValueError(
                    f"Molecule '{self.formula}': tag category must be a string."
                )
            if not isinstance(values, list) or not all(
                isinstance(value, str) for value in values
            ):
                raise ValueError(
                    f"Molecule '{self.formula}': tags['{category}'] "
                    "must be a list of strings."
                )

        # Validate ref roles.
        for role in self.refs:
            if role not in MOLECULE_REF_ROLES:
                raise ValueError(
                    f"Molecule '{self.formula}': "
                    f"unknown ref role '{role}'. "
                    f"Must be one of: {MOLECULE_REF_ROLES}"
                )

    @cached_property
    def _formula(self):
        """Parsed molmass Formula object used by formula-derived properties.

        This is cached because many properties below use the same parsed
        formula. The formula is parsed once per Molecule instance rather than
        reparsed every time mass, atom counts, charge, etc. are requested.

        The parser options are intentionally conservative: no group
        abbreviations, no oligo/polymer shorthand, no fractional stoichiometry,
        arithmetic allowed where molmass supports it, and no empty formulas.
        """
        return Formula(
            self.formula,
            parse_groups=False,
            parse_oligos=False,
            parse_fractions=False,
            parse_arithmetic=True,
            allow_empty=False,
        )

    @cached_property
    def _simple_formula_tokens(self):
        """Tokenize simple formula strings as a fallback for unsupported isotopes."""
        formula = self.formula.rstrip("+-")
        tokens = []
        position = 0

        while position < len(formula):
            match = SIMPLE_FORMULA_TOKEN.match(formula, position)
            if match is None:
                return None

            symbol = match.group("isotope_symbol") or match.group("symbol")
            isotope = match.group("isotope")
            count = int(match.group("count") or 1)
            tokens.append((symbol, isotope, count))
            position = match.end()

        return tokens

    def _simple_atom_counts(self, isotopic=False):
        """Return atom counts from the fallback tokenizer, if applicable."""
        if self._simple_formula_tokens is None:
            return None

        counts = {}
        for symbol, isotope, count in self._simple_formula_tokens:
            key = f"{isotope}{symbol}" if isotopic and isotope else symbol
            counts[key] = counts.get(key, 0) + count
        return counts

    def _fallback_formula_mass(self, average=False):
        """Compute mass when molmass lacks one explicitly requested isotope."""
        if self._simple_formula_tokens is None:
            raise FormulaError("unsupported fallback formula", self.formula, 0)

        total = 0.0
        for symbol, isotope, count in self._simple_formula_tokens:
            total += self._fallback_atom_mass(symbol, isotope, average) * count

        return total - self.charge * ELECTRON.mass

    def _fallback_atom_mass(self, symbol, isotope, average):
        """Return an atomic mass for fallback mass calculations."""
        if isotope:
            key = f"{isotope}{symbol}"
            if key in ISOTOPE_MASS_OVERRIDES:
                return ISOTOPE_MASS_OVERRIDES[key]
            return ELEMENTS[symbol].isotopes[int(isotope)].mass

        return ELEMENTS[symbol].mass if average else ELEMENTS[symbol].exactmass

    # === Computed properties: use formula via molmass, with explicit overrides ===

    @property
    def atom_counts(self):
        """Elemental composition with isotopic substitutions collapsed.

        Examples
        --------
        HCOOH -> {'H': 2, 'C': 1, 'O': 2}
        CH    -> {'C': 1, 'H': 1}

        Use isotope_counts if isotope labels need to be retained.
        """
        try:
            return {
                symbol: item.count
                for symbol, item in self._formula.composition(isotopic=False).items()
                if symbol != "e-"
            }
        except FormulaError:
            counts = self._simple_atom_counts(isotopic=False)
            if counts is None:
                raise
            return counts

    @property
    def atoms(self):
        """Legacy alias for atom_counts."""
        return self.atom_counts

    @property
    def isotope_counts(self):
        """Elemental composition retaining isotope labels where present."""
        try:
            return {
                symbol: item.count
                for symbol, item in self._formula.composition(isotopic=True).items()
                if symbol != "e-"
            }
        except FormulaError:
            counts = self._simple_atom_counts(isotopic=True)
            if counts is None:
                raise
            return counts

    @property
    def mass(self):
        """Exact/monoisotopic molecular mass in amu."""
        try:
            return self._formula.monoisotopic_mass
        except FormulaError:
            return self._fallback_formula_mass(average=False)

    @property
    def average_mass(self):
        """Average molecular mass from terrestrial isotopic abundances."""
        try:
            return self._formula.mass
        except FormulaError:
            return self._fallback_formula_mass(average=True)

    @property
    def nominal_mass(self):
        """Nominal integer mass of the formula."""
        try:
            return self._formula.nominal_mass
        except FormulaError:
            if self._simple_formula_tokens is None:
                raise

            total = 0
            for symbol, isotope, count in self._simple_formula_tokens:
                mass = int(isotope) if isotope else ELEMENTS[symbol].nominalmass
                total += mass * count
            return total

    @property
    def natoms(self):
        """Total number of nuclei in the formula."""
        try:
            return self._formula.atoms
        except FormulaError:
            counts = self._simple_atom_counts(isotopic=True)
            if counts is None:
                raise
            return sum(counts.values())

    @property
    def charge(self):
        """Formal charge of the molecule."""
        try:
            return self._formula.charge
        except FormulaError:
            signs = re.search(r"([+-]+)$", self.formula)
            if signs is None:
                return 0
            text = signs.group(1)
            return text.count("+") - text.count("-")

    @property
    def cation(self):
        """True if positively charged."""
        return self.charge > 0

    @property
    def anion(self):
        """True if negatively charged."""
        return self.charge < 0

    @property
    def neutral(self):
        """True if not charged."""
        return self.charge == 0

    @property
    def nelectrons(self):
        """Total number of electrons after correcting for formal charge.

        This is the sum of neutral atomic numbers minus the molecular charge.
        A cation has fewer electrons; an anion has more.
        """
        total = 0
        for symbol, count in self.atom_counts.items():
            total += ELEMENTS[symbol].number * count
        return total - self.charge

    @property
    def odd_electron(self):
        """True if the molecule has an odd number of electrons."""
        return self.nelectrons % 2 == 1

    @property
    def radical(self):
        """Whether the molecule should be treated as a radical.

        If radical_override is set, the curated value is returned. Otherwise,
        radical character is inferred from the odd-electron rule. 
        """
        if self.radical_override is not None:
            return self.radical_override
        return self.odd_electron

    @property
    def du(self):
        """Degree of unsaturation using the original astromol convention.

        This is only computed for molecules containing H/D, C, N, O, Cl, F,
        and S. It is a formula-based heuristic, not a structure assignment.
        Returns None outside the domain where the simple expression is useful.
        """
        inclusion = {"D", "C", "H", "N", "O", "Cl", "F", "S"}
        atoms = self.atoms
        if any(element not in inclusion for element in atoms):
            return None

        h = atoms.get("H", 0) + atoms.get("D", 0)
        c = atoms.get("C", 0)
        n = atoms.get("N", 0)
        cl = atoms.get("Cl", 0)
        f = atoms.get("F", 0)

        return 1 + 0.5 * (-h + 2 * c + n - cl - f)

    @property
    def maxdu(self):
        """Maximum degree of unsaturation for the same heavy-atom inventory.

        This preserves the original astromol convention and is only computed
        for the same limited element set as du.
        """
        inclusion = {"D", "C", "H", "N", "O", "Cl", "F", "S"}
        atoms = self.atoms
        if any(element not in inclusion for element in atoms):
            return None

        c = atoms.get("C", 0)
        n = atoms.get("N", 0)

        return 1 + 0.5 * (2 * c + n)

    # === Properties derived from spectroscopic constants ===

    @property
    def is_linear(self):
        """
        True if the molecule is linear.
        Determined by having only a B rotational constant (no A or C).
        Returns None if no rotational constants are available.
        """
        if self.rotcon is None:
            return None
        return (
            self.rotcon.A is None
            and self.rotcon.B is not None
            and self.rotcon.C is None
        )

    @property
    def kappa(self):
        """
        Ray's asymmetry parameter.
        Ranges from -1 (prolate) to +1 (oblate).
        Only computable if all three rotational constants are known.
        """
        if self.rotcon is None:
            return None
        A = self.rotcon.A
        B = self.rotcon.B
        C = self.rotcon.C
        if A is not None and B is not None and C is not None:
            if A != C:
                return (2 * B - A - C) / (A - C)
        return None

    def __repr__(self):
        """What you see when you print this object."""
        return f"{self.formula} ({self.name})"
    
@dataclass
class Detection:
    """A specific detection event of a molecule in an astronomical source."""

    # === Required fields ===
    molecule: str             # molecule label, e.g. "mol:CH"
                              # resolved to Molecule object during loading
    sources: list             # list of source nicks, e.g. ["SgrB2", "TMC1"]
                              # resolved to Source objects during loading
    telescopes: list          # list of telescope nicks, e.g. ["GBT", "IRAM30"]
                              # resolved to Telescope objects during loading
    wavelengths: list         # list from WAVELENGTHS, e.g. ["cm", "mm"]
    year: int                 # year of detection, e.g. 1937
    type: str                 # one of DETECTION_TYPES, e.g. "first"
    note: str = None

    # === Significance flag ===
    first: bool = False       # True if this is the first detection of this molecule
                              # in this context. Can be True for any type.

    # === Date precision — fill in what you know ===
    month: int = None         # 1-12, if known
    day: int = None           # 1-31, if known

    # === References ===
    # Keys are roles from DETECTION_REF_ROLES, values are lists of BibTeX keys.
    # e.g. {"observation": ["bib1", "bib2"]}
    # Resolved to lists of Ref objects during loading.
    refs: dict = None

    # === LaTeX generation ===
    latex_text: str = None    # sentence fragment for paper generation
                              # may contain {placeholder} syntax for dynamic resolution

    def __post_init__(self):
        """Validates types, wavelengths, ref roles, and sets defaults."""

        # Validate detection type
        if self.type not in DETECTION_TYPES:
            raise ValueError(
                f"Detection of '{self.molecule}': "
                f"unknown type '{self.type}'. "
                f"Must be one of: {DETECTION_TYPES}"
            )

        # Validate each wavelength in the list
        for wl in self.wavelengths:
            if wl not in WAVELENGTHS:
                raise ValueError(
                    f"Detection of '{self.molecule}': "
                    f"unknown wavelength '{wl}'. "
                    f"Must be one of: {WAVELENGTHS}"
                )

        # Default refs to empty dict
        if self.refs is None:
            self.refs = {}

        # Validate ref roles
        for role in self.refs:
            if role not in DETECTION_REF_ROLES:
                raise ValueError(
                    f"Detection of '{self.molecule}': "
                    f"unknown ref role '{role}'. "
                    f"Must be one of: {DETECTION_REF_ROLES}"
                )

    @property
    def sortdate(self):
        """
        Build a date for sorting/ordering purposes.
        Uses the 1st for any unknown month or day.
        """
        m = self.month if self.month is not None else 1
        d = self.day if self.day is not None else 1
        return date(self.year, m, d)

    def __repr__(self):
        """What you see when you print this object."""
        sources_str = ", ".join(
            s if isinstance(s, str) else s.nick
            for s in self.sources
        )
        first_str = " [FIRST]" if self.first else ""
        return f"{self.molecule} in {sources_str} ({self.year}, {self.type}{first_str})"
