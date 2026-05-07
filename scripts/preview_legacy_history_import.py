"""Preview import of record history from the legacy molecules.py file.

This writes staging artifacts only and does not modify production data.
"""

from __future__ import annotations

import ast
from collections import Counter
import json
from pathlib import Path
import re

from preview_latex_body_import import (
    TEX_PATH,
    extract_tex_subsections,
    molecule_indexes as latex_molecule_indexes,
    map_subsection_to_molecule,
)


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "astromol" / "data"
LEGACY_PATH = DATA / "molecules_legacy.py"
MOLECULES_PATH = DATA / "molecules.json"
PREVIEW_PATH = DATA / "molecules.history.preview.json"
IMPORT_PATH = DATA / "legacy_history_import.preview.json"
REPORT_PATH = DATA / "legacy_history_import_report.md"

LABEL_ALIASES = {
    "CNN1": "1-C10H7CN",
    "CNN2": "2-C10H7CN",
    "CAN1": "1-C12H7CN",
    "CAN5": "5-C12H7CN",
    "CNP1": "1-C16H9CN",
    "CNP2": "2-C16H9CN",
    "CNP4": "4-C16H9CN",
    "C9H7CN": "2-C9H7CN",
    "cC5H": "c-C5H",
    "cC5H4CCH2": "c-C5H4CCH2",
    "nCH3CH2CH2OH": "n-CH3CH2CH2OH",
    "iCH3CH2CH2OH": "i-CH3CH2CH2OH",
    "iC4H8": "i-C4H8",
}


def literal_or_none(node):
    try:
        return ast.literal_eval(node)
    except (TypeError, ValueError):
        return None


def keyword_value(call: ast.Call, name: str):
    for keyword in call.keywords:
        if keyword.arg == name:
            return literal_or_none(keyword.value)
    return None


def extract_legacy_entries() -> list[dict]:
    tree = ast.parse(LEGACY_PATH.read_text())
    entries = []
    for node in tree.body:
        if not isinstance(node, ast.Assign) or not isinstance(node.value, ast.Call):
            continue
        call = node.value
        if not isinstance(call.func, ast.Name) or call.func.id != "Molecule":
            continue
        target = node.targets[0]
        if not isinstance(target, ast.Name):
            continue
        change_log = keyword_value(call, "change_log") or {}
        census_version = keyword_value(call, "census_version")
        if not census_version and not change_log:
            continue
        entries.append(
            {
                "legacy_var": target.id,
                "legacy_line": node.lineno,
                "legacy_label": keyword_value(call, "label"),
                "astromol_name": keyword_value(call, "astromol_name"),
                "name": keyword_value(call, "name"),
                "formula": keyword_value(call, "formula"),
                "table_formula": keyword_value(call, "table_formula"),
                "census_version": census_version,
                "change_log": change_log,
            }
        )
    return entries


def unique_index(pairs: list[tuple[str, str]]) -> dict[str, str]:
    buckets = {}
    for key, value in pairs:
        if not key:
            continue
        buckets.setdefault(key, set()).add(value)
    return {key: next(iter(values)) for key, values in buckets.items() if len(values) == 1}


def molecule_indexes(molecules: list[dict]) -> dict[str, dict[str, str]]:
    labels = {row["label"].removeprefix("mol:"): row["label"] for row in molecules}
    return {
        "labels": labels,
        "table_formula": unique_index(
            [(row.get("table_formula"), row["label"]) for row in molecules]
        ),
        "formula": unique_index([(row.get("formula"), row["label"]) for row in molecules]),
        "name": unique_index(
            [((row.get("name") or "").lower(), row["label"]) for row in molecules]
        ),
    }


def map_legacy_entry(entry: dict, indexes: dict[str, dict[str, str]]) -> dict:
    candidates = []
    for method, raw in (
        ("legacy_label", entry.get("legacy_label")),
        ("astromol_name", entry.get("astromol_name")),
        ("legacy_var", entry.get("legacy_var")),
    ):
        if not raw:
            continue
        normalized = LABEL_ALIASES.get(raw, raw)
        if normalized in indexes["labels"]:
            candidates.append((method, indexes["labels"][normalized]))

    table_formula = entry.get("table_formula")
    if table_formula in indexes["table_formula"]:
        candidates.append(("table_formula", indexes["table_formula"][table_formula]))

    formula = entry.get("formula")
    if formula in indexes["formula"]:
        candidates.append(("formula", indexes["formula"][formula]))

    name = (entry.get("name") or "").lower()
    if name in indexes["name"]:
        candidates.append(("name", indexes["name"][name]))

    labels = list(dict.fromkeys(label for _, label in candidates))
    if len(labels) == 1:
        return {
            "status": "mapped",
            "label": labels[0],
            "mapping_methods": [method for method, label in candidates if label == labels[0]],
        }
    if len(labels) > 1:
        return {
            "status": "ambiguous",
            "label": None,
            "mapping_methods": [method for method, _ in candidates],
            "candidate_labels": labels,
        }
    return {"status": "unmapped", "label": None, "mapping_methods": []}


def census_from_version(version: str | None) -> str | None:
    if not version:
        return None
    match = re.match(r"^(\d{4})\.", str(version))
    return match.group(1) if match else None


def event_kind(summary: str, version: str) -> str:
    text = f"{summary} {version}".lower()
    if "initial" in text or "first entry" in text:
        return "added"
    if "correct" in text or "typo" in text or "fix" in text:
        return "corrected"
    return "updated"


def event_fields(summary: str) -> list[str]:
    text = summary.lower()
    fields = []
    checks = [
        ("tag", "tags"),
        ("functional group", "tags"),
        ("detection", "detections"),
        ("detected", "detections"),
        ("isotopologue", "detections"),
        ("ppd", "detections"),
        ("ice", "detections"),
        ("exgal", "detections"),
        ("exoplanet", "detections"),
        ("reference", "refs"),
        ("ref", "refs"),
        ("lab", "refs"),
        ("name", "name"),
        ("formula", "formula"),
        ("note", "note"),
    ]
    for needle, field in checks:
        if needle in text and field not in fields:
            fields.append(field)
    return fields


def sorted_change_items(change_log: dict, census_version: str | None) -> list[tuple[str, str]]:
    items = []
    if change_log:
        for version, summary in change_log.items():
            items.append((str(version), str(summary)))
    elif census_version:
        items.append((str(census_version), "Initial entry"))

    def key(item):
        version = item[0]
        parts = []
        for part in version.split("."):
            try:
                parts.append(int(part))
            except ValueError:
                parts.append(0)
        return parts

    return sorted(items, key=key)


def labels_in_2021_tex(molecules: list[dict]) -> set[str]:
    """Return current molecule labels represented in the 2021 manuscript TeX."""
    indexes = latex_molecule_indexes(molecules)
    labels = set()
    for subsection in extract_tex_subsections(TEX_PATH.read_text()):
        mapping = map_subsection_to_molecule(subsection, indexes)
        if mapping["status"] == "mapped":
            labels.add(mapping["label"])
    return labels


def print_census_for_legacy_version(version: str, baseline_census: str) -> str | None:
    """Map legacy database-version strings to semantic print-census releases."""
    version_census = census_from_version(version)
    if baseline_census == "2026" and version_census in {"2018", "2021"}:
        return "2026"
    return version_census


def build_history(entry: dict, baseline_census: str) -> tuple[dict, list[str]]:
    events = []
    warnings = []
    for version, summary in sorted_change_items(
        entry.get("change_log") or {}, entry.get("census_version")
    ):
        census = print_census_for_legacy_version(version, baseline_census)
        if census not in {"2018", "2021", "2026"}:
            warnings.append(f"unexpected census prefix `{census}` in `{version}`")
        events.append(
            {
                "date": None,
                "census": census,
                "kind": event_kind(summary, version),
                "fields": event_fields(summary),
                "summary": summary,
                "legacy_version": version,
            }
        )

    introduced_event = next((event for event in events if event["kind"] == "added"), None)
    if introduced_event is None and events:
        introduced_event = events[0]
    introduced = {}
    if introduced_event:
        introduced = {
            "date": introduced_event["date"],
            "census": introduced_event["census"],
            "legacy_version": introduced_event["legacy_version"],
        }

    return (
        {
            "introduced": introduced,
            "last_modified": None,
            "last_reviewed": None,
            "events": events,
        },
        warnings,
    )


def make_report(rows: list[dict]) -> str:
    counts = Counter(row["status"] for row in rows)
    event_counts = Counter(
        event["kind"]
        for row in rows
        if row["status"] == "mapped"
        for event in row["history"]["events"]
    )
    warning_rows = [row for row in rows if row.get("warnings")]
    problem_rows = [row for row in rows if row["status"] != "mapped"]
    new_2026_rows = [
        row for row in rows
        if row["status"] == "mapped"
        and row["history"]["introduced"].get("census") == "2026"
    ]

    lines = []
    lines.append("# Legacy History Import Preview")
    lines.append("")
    lines.append("Generated by `scripts/preview_legacy_history_import.py`.")
    lines.append("")
    lines.append("## Summary")
    lines.append("")
    lines.append(f"- Legacy molecule histories extracted: {len(rows)}")
    lines.append(f"- Mapped histories: {counts.get('mapped', 0)}")
    lines.append(f"- Ambiguous histories: {counts.get('ambiguous', 0)}")
    lines.append(f"- Unmapped histories: {counts.get('unmapped', 0)}")
    lines.append(f"- `added` events: {event_counts.get('added', 0)}")
    lines.append(f"- `updated` events: {event_counts.get('updated', 0)}")
    lines.append(f"- `corrected` events: {event_counts.get('corrected', 0)}")
    lines.append(f"- Histories introduced in the 2026 print census: {len(new_2026_rows)}")
    lines.append("")
    lines.append("Preview artifacts:")
    lines.append(f"- `{PREVIEW_PATH.relative_to(ROOT)}`: molecules with staged `history` values.")
    lines.append(f"- `{IMPORT_PATH.relative_to(ROOT)}`: per-legacy-entry mapping details.")
    lines.append("")

    lines.append("## Mapping Issues")
    lines.append("")
    if not problem_rows:
        lines.append("No mapping issues found.")
    else:
        for row in problem_rows:
            lines.append(
                f"- line {row['legacy_line']} `{row['legacy_var']}` "
                f"formula `{row.get('formula')}` name `{row.get('name')}` -> {row['status']}"
            )
            if row.get("candidate_labels"):
                lines.append(f"  - candidates: {', '.join(row['candidate_labels'])}")
    lines.append("")

    lines.append("## Warnings")
    lines.append("")
    if not warning_rows:
        lines.append("No warnings.")
    else:
        for row in warning_rows:
            lines.append(f"- line {row['legacy_line']} `{row['legacy_var']}` -> `{row.get('label')}`")
            for warning in row["warnings"]:
                lines.append(f"  - {warning}")
    lines.append("")

    lines.append("## Introduced In 2026 Print Census")
    lines.append("")
    if not new_2026_rows:
        lines.append("No legacy histories were classified as first appearing in the 2026 print census.")
    else:
        lines.append(
            "These legacy records were absent from `2021_census_arxiv.tex`, so their "
            "`history.introduced.census` is staged as `2026` even when the preserved "
            "`legacy_version` begins with `2021`."
        )
        lines.append("")
        for row in new_2026_rows:
            intro = row["history"]["introduced"]
            lines.append(
                f"- line {row['legacy_line']} `{row['legacy_var']}` -> `{row['label']}` "
                f"introduced `{intro.get('census')}` from legacy `{intro.get('legacy_version')}`"
            )
    lines.append("")

    lines.append("## Mapped History Events")
    lines.append("")
    for row in [row for row in rows if row["status"] == "mapped"]:
        lines.append(
            f"### line {row['legacy_line']} `{row['legacy_var']}` -> `{row['label']}`"
        )
        lines.append("")
        lines.append(f"- mapping methods: `{', '.join(row['mapping_methods'])}`")
        lines.append(f"- introduced: `{row['history']['introduced']}`")
        lines.append("")
        for event in row["history"]["events"]:
            fields = ", ".join(event["fields"]) if event["fields"] else ""
            lines.append(
                f"- `{event['legacy_version']}` `{event['kind']}` "
                f"census `{event['census']}` fields `{fields}`: {event['summary']}"
            )
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def main() -> None:
    molecules = json.loads(MOLECULES_PATH.read_text())
    indexes = molecule_indexes(molecules)
    preview = json.loads(json.dumps(molecules))
    preview_by_label = {row["label"]: row for row in preview}
    baseline_2021_labels = labels_in_2021_tex(molecules)

    rows = []
    for entry in extract_legacy_entries():
        mapping = map_legacy_entry(entry, indexes)
        baseline_census = (
            "2021" if mapping["label"] in baseline_2021_labels else "2026"
        )
        history, warnings = build_history(entry, baseline_census)
        row = {
            **entry,
            **mapping,
            "baseline_census": baseline_census,
            "history": history,
            "warnings": warnings,
        }
        if mapping["status"] == "mapped":
            preview_by_label[mapping["label"]]["history"] = history
        rows.append(row)

    PREVIEW_PATH.write_text(json.dumps(preview, indent=2) + "\n")
    IMPORT_PATH.write_text(json.dumps(rows, indent=2) + "\n")
    REPORT_PATH.write_text(make_report(rows))

    print(f"Wrote {PREVIEW_PATH.relative_to(ROOT)}")
    print(f"Wrote {IMPORT_PATH.relative_to(ROOT)}")
    print(f"Wrote {REPORT_PATH.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
