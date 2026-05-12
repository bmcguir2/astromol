"""Validation helpers for astromol production data."""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import date
import numbers

from .database import Database
from .models import DETECTION_RELATION_FIELDS


RELATION_RECIPROCALS = {
    "confirms": "confirmed_by",
    "confirmed_by": "confirms",
    "disputes": "disputed_by",
    "disputed_by": "disputes",
    "supersedes": "superseded_by",
    "superseded_by": "supersedes",
}


@dataclass(frozen=True)
class ValidationIssue:
    """One validation issue found in the production data."""

    severity: str
    code: str
    record: str
    message: str


@dataclass
class ValidationReport:
    """Collection of validation issues."""

    issues: list[ValidationIssue] = field(default_factory=list)

    @property
    def errors(self) -> list[ValidationIssue]:
        return [issue for issue in self.issues if issue.severity == "error"]

    @property
    def warnings(self) -> list[ValidationIssue]:
        return [issue for issue in self.issues if issue.severity == "warning"]

    @property
    def ok(self) -> bool:
        return not self.errors

    def add(self, severity: str, code: str, record: str, message: str) -> None:
        self.issues.append(
            ValidationIssue(
                severity=severity,
                code=code,
                record=record,
                message=message,
            )
        )

    def raise_for_errors(self) -> None:
        if self.errors:
            details = "\n".join(format_issue(issue) for issue in self.errors)
            raise ValueError(f"astromol validation failed:\n{details}")


def validate_database(db: Database | None = None) -> ValidationReport:
    """Validate production records and return a report.

    ``Database()`` already validates parsing, required dataclass values, duplicate
    primary keys, reference resolution, and relationship targets. This layer adds
    semantic consistency checks that are useful before releases and data-update
    commits.
    """
    if db is None:
        db = Database()

    report = ValidationReport()
    _validate_molecule_labels(db, report)
    _validate_detection_ids(db, report)
    _validate_detection_relationship_reciprocals(db, report)
    _validate_history(db, report)
    _validate_spectroscopy_values(db, report)
    _validate_first_detection_flags(db, report)
    return report


def _validate_molecule_labels(db: Database, report: ValidationReport) -> None:
    for molecule in db.molecules.values():
        if not molecule.label.startswith("mol:"):
            report.add(
                "error",
                "molecule-label",
                molecule.label,
                "Molecule labels must start with 'mol:'.",
            )


def _validate_detection_ids(db: Database, report: ValidationReport) -> None:
    for detection in db.detections:
        if not detection.id.startswith("det:"):
            report.add(
                "error",
                "detection-id",
                detection.id,
                "Detection IDs must start with 'det:'.",
            )


def _validate_detection_relationship_reciprocals(
    db: Database,
    report: ValidationReport,
) -> None:
    for detection in db.detections:
        for field_name in DETECTION_RELATION_FIELDS:
            reciprocal_name = RELATION_RECIPROCALS[field_name]
            for target_id in getattr(detection, field_name):
                target = db.detections_by_id[target_id]
                if detection.id not in getattr(target, reciprocal_name):
                    report.add(
                        "error",
                        "relationship-reciprocal",
                        detection.id,
                        (
                            f"{field_name} points to {target_id}, but "
                            f"{target_id}.{reciprocal_name} does not include "
                            f"{detection.id}."
                        ),
                    )


def _validate_history(db: Database, report: ValidationReport) -> None:
    for kind, records in (
        ("molecule", db.molecules.values()),
        ("detection", db.detections),
        ("source", db.sources.values()),
        ("telescope", db.telescopes.values()),
    ):
        for record in records:
            record_id = getattr(record, "label", None) or getattr(record, "id", None) or getattr(record, "nick", "")
            history = getattr(record, "history", None)
            if history is None:
                if kind in {"molecule", "detection"}:
                    report.add(
                        "warning",
                        "history-missing",
                        record_id,
                        f"{kind} record has no history metadata.",
                    )
                continue

            introduced = _census_int(history.introduced_census)
            accepted = _census_int(history.accepted_census)
            if accepted is not None and introduced is not None and accepted < introduced:
                report.add(
                    "error",
                    "history-order",
                    record_id,
                    (
                        "history.accepted.census is earlier than "
                        "history.introduced.census."
                    ),
                )
            if accepted is not None and history.accepted_context is None:
                report.add(
                    "warning",
                    "history-accepted-context",
                    record_id,
                    "history.accepted is missing a context value.",
                )
            if history.last_modified is not None:
                try:
                    date.fromisoformat(history.last_modified)
                except ValueError:
                    report.add(
                        "error",
                        "history-date",
                        record_id,
                        "history.last_modified must use ISO format YYYY-MM-DD.",
                    )


def _validate_spectroscopy_values(db: Database, report: ValidationReport) -> None:
    for molecule in db.molecules.values():
        if molecule.rotcon is not None:
            for component in ("A", "B", "C"):
                value = getattr(molecule.rotcon, component)
                if value is not None and not _is_number(value):
                    report.add(
                        "error",
                        "rotcon-nonnumeric",
                        molecule.label,
                        f"rotcon.{component} must be numeric or null; found {value!r}.",
                    )
        if molecule.dipole is not None:
            for component in ("a", "b", "c"):
                value = getattr(molecule.dipole, component)
                if value is None or _is_number(value):
                    continue
                severity = "warning" if value == "*" else "error"
                report.add(
                    severity,
                    "dipole-nonnumeric",
                    molecule.label,
                    f"dipole.{component} should be numeric or null; found {value!r}.",
                )


def _validate_first_detection_flags(
    db: Database,
    report: ValidationReport,
) -> None:
    groups = {}
    for detection in db.detections:
        if detection.status != "secure":
            continue
        if detection.molecule.isotopologue_of is not None:
            continue
        key = (detection.molecule.label, detection.type)
        groups.setdefault(key, []).append(detection)

    for (molecule_label, detection_type), detections in groups.items():
        expected = min(
            detections,
            key=lambda detection: (
                detection.sortdate,
                detection.id,
            ),
        )
        first_flags = [detection for detection in detections if detection.first]
        if not first_flags:
            continue
        if expected not in first_flags:
            report.add(
                "error",
                "first-flag-wrong-record",
                molecule_label,
                (
                    f"{expected.id} is the earliest secure {detection_type} "
                    "detection, but it is not marked first."
                ),
            )
        for detection in first_flags:
            if detection is not expected:
                report.add(
                    "error",
                    "first-flag-extra",
                    detection.id,
                    (
                        f"Detection is marked first, but {expected.id} is the "
                        f"earliest secure {detection_type} detection for "
                        f"{molecule_label}."
                    ),
                )


def _is_number(value) -> bool:
    return isinstance(value, numbers.Real) and not isinstance(value, bool)


def _census_int(value: str | int | None) -> int | None:
    if value is None:
        return None
    return int(value)


def format_issue(issue: ValidationIssue) -> str:
    """Return one issue formatted for command-line output."""
    return f"[{issue.severity}] {issue.code} {issue.record}: {issue.message}"


def format_report(report: ValidationReport) -> str:
    """Return a human-readable validation report."""
    lines = [
        "# astromol validation report",
        "",
        f"- errors: {len(report.errors)}",
        f"- warnings: {len(report.warnings)}",
    ]
    if report.issues:
        lines.append("")
        for issue in report.issues:
            lines.append(format_issue(issue))
    return "\n".join(lines)


def main() -> int:
    report = validate_database()
    print(format_report(report))
    return 0 if report.ok else 1


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
