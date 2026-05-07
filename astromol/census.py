"""Census-scoped views over an :class:`astromol.database.Database`.

The view layer centralizes record selection for manuscript tables, figures,
slides, and statistics. Output code should consume these views instead of
reimplementing census/history filtering locally.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from typing import Iterable

from .models import DETECTION_TYPES, Detection, Molecule, RecordHistory


CURRENT_SCOPE = "current"
CENSUS_SCOPE = "census"


def _census_int(census: str | int | None) -> int | None:
    """Normalize census labels such as ``"2021"`` to integers."""
    if census is None:
        return None
    return int(census)


def _history_census(history: RecordHistory | None, field: str) -> int | None:
    """Return an integer census value from a RecordHistory property."""
    if history is None:
        return None
    value = getattr(history, field)
    return _census_int(value)


@dataclass(frozen=True)
class CensusView:
    """A filtered view of astromol records for a census boundary or live data.

    Use :meth:`for_census` for frozen/historical census outputs and
    :meth:`current` for the live database. In census mode, accepted records are
    selected by ``history.accepted.census <= census``. In current mode,
    accepted records are selected regardless of census label. Context queries
    exclude isotopologues by default; pass ``include_isotopologues=True`` for
    isotope-expanded views.
    """

    db: object
    scope: str = CENSUS_SCOPE
    census: str | None = None

    def __post_init__(self):
        if self.scope not in {CENSUS_SCOPE, CURRENT_SCOPE}:
            raise ValueError(
                f"Unknown CensusView scope '{self.scope}'. "
                f"Use '{CENSUS_SCOPE}' or '{CURRENT_SCOPE}'."
            )
        if self.scope == CENSUS_SCOPE and self.census is None:
            raise ValueError("CensusView census mode requires a census value.")
        if self.census is not None:
            object.__setattr__(self, "census", str(self.census))

    @classmethod
    def for_census(cls, db, census: str | int) -> "CensusView":
        """Create a view frozen at a census boundary."""
        return cls(db=db, scope=CENSUS_SCOPE, census=str(census))

    @classmethod
    def current(cls, db) -> "CensusView":
        """Create a live view of all currently accepted records."""
        return cls(db=db, scope=CURRENT_SCOPE, census=None)

    @property
    def is_current(self) -> bool:
        """Whether this view represents live database contents."""
        return self.scope == CURRENT_SCOPE

    @property
    def census_year(self) -> int | None:
        """The integer census boundary, or ``None`` for current views."""
        return _census_int(self.census)

    def accepted_record(self, record) -> bool:
        """Return whether a molecule or detection is accepted in this view."""
        history = getattr(record, "history", None)
        if history is None or history.accepted is None:
            return False
        if self.is_current:
            return True
        accepted = _history_census(history, "accepted_census")
        return accepted is not None and accepted <= self.census_year

    def introduced_record(self, record) -> bool:
        """Return whether a record had entered tracking by this view."""
        history = getattr(record, "history", None)
        if history is None:
            return False
        if self.is_current:
            return bool(history.introduced)
        introduced = _history_census(history, "introduced_census")
        return introduced is not None and introduced <= self.census_year

    @staticmethod
    def is_secure(detection: Detection) -> bool:
        """Return whether a detection is treated as secure."""
        return detection.status == "secure"

    def detection_in_scope(
        self,
        detection: Detection,
        *,
        include_tentative: bool = False,
        include_disputed: bool = False,
    ) -> bool:
        """Return whether a detection belongs in this view.

        Secure detections require accepted history. Tentative and disputed
        detections are included only when requested, and then by introduced
        history rather than accepted history.
        """
        if self.is_secure(detection):
            return self.accepted_record(detection)
        if detection.status == "tentative" and include_tentative:
            return self.introduced_record(detection)
        if detection.status == "disputed" and include_disputed:
            return self.introduced_record(detection)
        return False

    def detections(
        self,
        detection_type: str | None = None,
        *,
        include_tentative: bool = False,
        include_disputed: bool = False,
        include_isotopologues: bool = False,
    ) -> list[Detection]:
        """Return detections in this view, optionally restricted by type.

        Isotopologue detections are excluded by default. Set
        ``include_isotopologues=True`` when the expanded isotopologue inventory
        is wanted.
        """
        if detection_type is not None and detection_type not in DETECTION_TYPES:
            raise ValueError(
                f"Unknown detection type '{detection_type}'. "
                f"Must be one of: {DETECTION_TYPES}"
            )

        detections = [
            detection
            for detection in self.db.detections
            if (detection_type is None or detection.type == detection_type)
            and (
                include_isotopologues
                or detection.molecule.isotopologue_of is None
            )
            and self.detection_in_scope(
                detection,
                include_tentative=include_tentative,
                include_disputed=include_disputed,
            )
        ]
        return sorted(
            detections,
            key=lambda detection: (
                detection.sortdate,
                detection.molecule.label,
                detection.id,
            ),
        )

    def context_detections(
        self,
        detection_type: str,
        *,
        include_tentative: bool = False,
        include_disputed: bool = False,
        include_isotopologues: bool = False,
    ) -> list[Detection]:
        """Return detections for one context such as ``"ISM/CSM"``."""
        return self.detections(
            detection_type,
            include_tentative=include_tentative,
            include_disputed=include_disputed,
            include_isotopologues=include_isotopologues,
        )

    def context_molecules(
        self,
        detection_type: str,
        *,
        include_tentative: bool = False,
        include_disputed: bool = False,
        include_isotopologues: bool = False,
    ) -> list[Molecule]:
        """Return unique molecules detected in one context.

        Isotopologues are excluded by default and can be included with
        ``include_isotopologues=True``.
        """
        molecules = {
            detection.molecule.label: detection.molecule
            for detection in self.context_detections(
                detection_type,
                include_tentative=include_tentative,
                include_disputed=include_disputed,
                include_isotopologues=include_isotopologues,
            )
        }
        return sorted(
            molecules.values(),
            key=lambda molecule: (
                molecule.natoms,
                molecule.label,
            ),
        )

    def accepted_molecules(
        self,
        *,
        include_isotopologues: bool = False,
    ) -> list[Molecule]:
        """Return molecule records accepted in this view by molecule history.

        Isotopologues are excluded by default and can be included with
        ``include_isotopologues=True``.
        """
        return sorted(
            (
                molecule
                for molecule in self.db.molecules.values()
                if self.accepted_record(molecule)
                and (
                    include_isotopologues
                    or molecule.isotopologue_of is None
                )
            ),
            key=lambda molecule: (
                molecule.natoms,
                molecule.label,
            ),
        )

    def ism_detections(self, **kwargs) -> list[Detection]:
        """Return ISM/CSM detections in this view."""
        return self.context_detections("ISM/CSM", **kwargs)

    def ism_molecules(self, **kwargs) -> list[Molecule]:
        """Return ISM/CSM molecules in this view."""
        return self.context_molecules("ISM/CSM", **kwargs)

    def ppd_detections(self, **kwargs) -> list[Detection]:
        """Return protoplanetary-disk detections in this view."""
        return self.context_detections("ppd", **kwargs)

    def ppd_molecules(self, **kwargs) -> list[Molecule]:
        """Return protoplanetary-disk molecules in this view."""
        return self.context_molecules("ppd", **kwargs)

    def ice_detections(self, **kwargs) -> list[Detection]:
        """Return ice detections in this view."""
        return self.context_detections("ice", **kwargs)

    def ice_molecules(self, **kwargs) -> list[Molecule]:
        """Return ice molecules in this view."""
        return self.context_molecules("ice", **kwargs)

    def exgal_detections(self, **kwargs) -> list[Detection]:
        """Return extragalactic detections in this view."""
        return self.context_detections("exgal", **kwargs)

    def exgal_molecules(self, **kwargs) -> list[Molecule]:
        """Return extragalactic molecules in this view."""
        return self.context_molecules("exgal", **kwargs)

    def exoplanet_detections(self, **kwargs) -> list[Detection]:
        """Return exoplanet-atmosphere detections in this view."""
        return self.context_detections("exo", **kwargs)

    def exoplanet_molecules(self, **kwargs) -> list[Molecule]:
        """Return exoplanet-atmosphere molecules in this view."""
        return self.context_molecules("exo", **kwargs)

    def source_counts(
        self,
        detection_type: str = "ISM/CSM",
        *,
        include_tentative: bool = False,
        include_disputed: bool = False,
        include_isotopologues: bool = False,
        key: str = "nick",
        group_diffuse_cloud: bool = False,
        diffuse_cloud_label: str = "DiffuseCloud",
    ) -> Counter:
        """Count source contributions for detections in this view.

        ``key`` may be ``"nick"``, ``"name"``, or ``"latex_name"``.
        Set ``group_diffuse_cloud=True`` to consolidate all diffuse-cloud/LOS
        sources under one label for historical manuscript table reproduction.
        """
        if key not in {"nick", "name", "latex_name"}:
            raise ValueError("source_counts key must be 'nick', 'name', or 'latex_name'.")

        counts = Counter()
        for detection in self.context_detections(
            detection_type,
            include_tentative=include_tentative,
            include_disputed=include_disputed,
            include_isotopologues=include_isotopologues,
        ):
            for source in detection.sources:
                if group_diffuse_cloud and source.type == "Diffuse Cloud":
                    counts[diffuse_cloud_label] += 1
                elif key == "latex_name":
                    counts[source.latex_name or source.name] += 1
                else:
                    counts[getattr(source, key)] += 1
        return counts

    def facility_counts(
        self,
        detection_type: str = "ISM/CSM",
        *,
        include_tentative: bool = False,
        include_disputed: bool = False,
        include_isotopologues: bool = False,
        key: str = "nick",
    ) -> Counter:
        """Count telescope/facility contributions for detections in this view."""
        if key not in {"nick", "shortname", "name", "latex_name"}:
            raise ValueError(
                "facility_counts key must be 'nick', 'shortname', 'name', or 'latex_name'."
            )

        counts = Counter()
        for detection in self.context_detections(
            detection_type,
            include_tentative=include_tentative,
            include_disputed=include_disputed,
            include_isotopologues=include_isotopologues,
        ):
            for telescope in detection.telescopes:
                if key == "latex_name":
                    counts[telescope.latex_name or telescope.shortname or telescope.name] += 1
                else:
                    counts[getattr(telescope, key)] += 1
        return counts

    def molecule_labels(self, molecules: Iterable[Molecule]) -> set[str]:
        """Return molecule labels from a molecule iterable."""
        return {molecule.label for molecule in molecules}
