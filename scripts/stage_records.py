"""Stage compact YAML curation records into preview JSON.

The production data remains JSON. This script lets curators write small YAML
records, validates them against the current model, and writes preview artifacts.
Use --apply only after reviewing the generated report.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from datetime import date
import json
from pathlib import Path
import re
import sys

try:
    import yaml
except ImportError as exc:  # pragma: no cover - environment guard
    raise SystemExit("stage_records.py requires PyYAML: pip install pyyaml") from exc


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "astromol" / "data"
CURRENT_CENSUS = "2026"

sys.path.insert(0, str(ROOT))

from astromol.models import (  # noqa: E402
    Detection,
    Molecule,
    Source,
    Telescope,
    DETECTION_RELATION_FIELDS,
    DETECTION_REF_ROLES,
    MOLECULE_REF_ROLES,
)


KIND_TO_FILE = {
    "molecule": "molecules.json",
    "detection": "detections.json",
    "source": "sources.json",
    "telescope": "telescopes.json",
}

FIELDS = {
    "molecule": [
        "name",
        "formula",
        "table_formula",
        "label",
        "note",
        "iupac_name",
        "selfies",
        "synonyms",
        "smiles",
        "canonical_smiles",
        "inchi",
        "inchikey",
        "radical_override",
        "fullerene",
        "pah",
        "n_rings",
        "cyclic",
        "tags",
        "rotcon",
        "dipole",
        "refs",
        "isotopologue_of",
        "latex_section_override",
        "latex_body",
        "history",
    ],
    "detection": [
        "id",
        "note",
        "molecule",
        "sources",
        "telescopes",
        "wavelengths",
        "year",
        "type",
        "status",
        "status_note",
        "first",
        "month",
        "day",
        "refs",
        "confirms",
        "confirmed_by",
        "disputes",
        "disputed_by",
        "supersedes",
        "superseded_by",
        "latex_text",
        "history",
    ],
    "source": [
        "name",
        "nick",
        "type",
        "ra",
        "dec",
        "simbad_url",
        "latex_name",
        "note",
        "history",
    ],
    "telescope": [
        "name",
        "nick",
        "shortname",
        "type",
        "wavelength",
        "diameter",
        "latitude",
        "longitude",
        "built",
        "decommissioned",
        "note",
        "latex_name",
        "history",
    ],
}

DEFAULTS = {
    "molecule": {
        "table_formula": None,
        "note": None,
        "iupac_name": None,
        "selfies": None,
        "synonyms": [],
        "smiles": None,
        "canonical_smiles": None,
        "inchi": None,
        "inchikey": None,
        "radical_override": None,
        "fullerene": False,
        "pah": False,
        "n_rings": 0,
        "cyclic": False,
        "tags": {},
        "rotcon": None,
        "dipole": None,
        "refs": {"lab": [], "computation": []},
        "isotopologue_of": None,
        "latex_section_override": None,
        "latex_body": None,
        "history": None,
    },
    "detection": {
        "note": None,
        "sources": [],
        "telescopes": [],
        "wavelengths": [],
        "status": "secure",
        "status_note": None,
        "first": False,
        "month": None,
        "day": None,
        "refs": {},
        "confirms": [],
        "confirmed_by": [],
        "disputes": [],
        "disputed_by": [],
        "supersedes": [],
        "superseded_by": [],
        "latex_text": None,
        "history": None,
    },
    "source": {
        "ra": None,
        "dec": None,
        "simbad_url": None,
        "latex_name": None,
        "note": None,
        "history": None,
    },
    "telescope": {
        "diameter": None,
        "latitude": None,
        "longitude": None,
        "built": None,
        "decommissioned": None,
        "note": None,
        "latex_name": None,
        "history": None,
    },
}

REQUIRED = {
    "molecule": ["name", "formula", "label"],
    "detection": ["id", "molecule", "year", "type"],
    "source": ["name", "nick", "type"],
    "telescope": ["name", "nick", "shortname", "type", "wavelength"],
}

DATACLASS = {
    "molecule": Molecule,
    "detection": Detection,
    "source": Source,
    "telescope": Telescope,
}

REF_ROLE_FIELDS = {
    "molecule": MOLECULE_REF_ROLES,
    "detection": DETECTION_REF_ROLES,
}

MISSING = object()


def load_json(name: str) -> list[dict]:
    return json.loads((DATA / name).read_text())


def write_json(path: Path, data: list[dict]) -> None:
    path.write_text(json.dumps(data, indent=2) + "\n")


def load_reference_keys() -> set[str]:
    text = (DATA / "references.bib").read_text()
    return set(re.findall(r"@\w+\{([^,\n]+)", text))


def load_yaml_records(paths: list[Path]) -> list[dict]:
    records = []
    for path in paths:
        loaded = yaml.safe_load(path.read_text())
        if loaded is None:
            continue
        if isinstance(loaded, dict) and "records" in loaded:
            loaded = loaded["records"]
        if not isinstance(loaded, list):
            raise ValueError(f"{path}: expected a YAML list or a dict with records: [...]")
        for index, record in enumerate(loaded, start=1):
            if not isinstance(record, dict):
                raise ValueError(f"{path}: record {index} is not a mapping")
            record = dict(record)
            record["_staging_file"] = str(path)
            record["_staging_index"] = index
            records.append(record)
    return records


def prune_empty_template_value(value):
    """Return MISSING for untouched optional template values."""
    if value is None:
        return MISSING
    if isinstance(value, str) and value == "":
        return MISSING
    if isinstance(value, list):
        kept = []
        for item in value:
            pruned = prune_empty_template_value(item)
            if pruned is not MISSING:
                kept.append(pruned)
        return kept if kept else MISSING
    if isinstance(value, dict):
        kept = {}
        for key, item in value.items():
            pruned = prune_empty_template_value(item)
            if pruned is not MISSING:
                kept[key] = pruned
        return kept if kept else MISSING
    return value


def clean_payload_value(key: str, value):
    """Prune blank template values without dropping meaningful false/zero values."""
    if key == "refs" and isinstance(value, dict):
        refs = {}
        for role, keys in value.items():
            if keys:
                refs[role] = keys
        return refs if refs else MISSING

    if key == "history" and isinstance(value, dict):
        history = {}
        for history_key, item in value.items():
            if history_key == "accepted" and item is None:
                history[history_key] = None
                continue
            pruned = prune_empty_template_value(item)
            if pruned is not MISSING:
                history[history_key] = pruned
        return history if history else MISSING

    if key in {"rotcon", "dipole", "history", "tags"}:
        return prune_empty_template_value(value)

    return prune_empty_template_value(value)


def default_history(kind: str, run_date: str) -> dict:
    introduced = {
        "date": run_date,
    }
    history = {
        "introduced": introduced,
        "last_modified": run_date,
        "events": [
            {
                "kind": "added",
                "summary": f"Added {kind} record.",
                "date": run_date,
                "fields": [],
            }
        ],
    }
    if kind == "molecule":
        introduced["context"] = "confirmed"
        history["accepted"] = {
            "date": run_date,
            "census": CURRENT_CENSUS,
            "context": "confirmed",
        }
    return history


def normalize_history(kind: str, history: dict | None, run_date: str) -> dict:
    if history is None:
        return default_history(kind, run_date)

    history = dict(history)
    history.pop("last_reviewed", None)
    introduced = dict(history.get("introduced") or {})
    if not introduced.get("date"):
        introduced["date"] = run_date
    if kind == "molecule" and not introduced.get("context"):
        introduced["context"] = (
            "tentative" if "accepted" in history and history["accepted"] is None
            else "confirmed"
        )
    history["introduced"] = introduced

    if kind == "molecule":
        if "accepted" not in history:
            history["accepted"] = default_history(kind, run_date)["accepted"]
        elif history["accepted"] is not None:
            accepted = dict(history["accepted"] or {})
            if not accepted.get("date"):
                accepted["date"] = run_date
            if not accepted.get("census"):
                accepted["census"] = CURRENT_CENSUS
            if not accepted.get("context"):
                accepted["context"] = "confirmed"
            history["accepted"] = accepted
    else:
        history.pop("accepted", None)

    if not history.get("last_modified"):
        history["last_modified"] = run_date

    events = []
    for event in history.get("events") or []:
        event = dict(event)
        if not event.get("date"):
            event["date"] = run_date
        events.append(event)
    if not events:
        events = default_history(kind, run_date)["events"]
    history["events"] = events
    return history


def merge_defaults(kind: str, payload: dict) -> dict:
    merged = json.loads(json.dumps(DEFAULTS.get(kind, {})))
    for key, value in payload.items():
        if key == "kind" or key.startswith("_"):
            continue
        value = clean_payload_value(key, value)
        if value is MISSING:
            continue
        if key == "refs" and isinstance(value, dict):
            refs = dict(merged.get("refs") or {})
            refs.update(value)
            merged["refs"] = refs
        else:
            merged[key] = value

    if kind == "molecule" and merged.get("table_formula") is None:
        merged["table_formula"] = merged.get("formula")

    return {field: merged.get(field) for field in FIELDS[kind]}


def staging_meta(record: dict) -> dict:
    return {
        key: value
        for key, value in record.items()
        if key.startswith("_") and key not in {"_staging_file", "_staging_index"}
    }


def ref_keys_from(payload: dict) -> list[str]:
    refs = payload.get("refs") or {}
    return [key for keys in refs.values() for key in keys]


def validate_record(
    kind: str,
    payload: dict,
    existing: dict[str, set[str]],
    staged: dict[str, set[str]],
    ref_keys: set[str],
) -> list[str]:
    errors = []

    for field in REQUIRED[kind]:
        value = payload.get(field)
        if value is None or value == "" or value == [] or value == {}:
            errors.append(f"missing required field `{field}`")

    unknown = [
        key
        for key in payload
        if key not in FIELDS[kind]
    ]
    if unknown:
        errors.append(f"unknown output fields: {', '.join(unknown)}")

    try:
        DATACLASS[kind](**payload)
    except Exception as exc:  # noqa: BLE001 - report validation detail
        errors.append(str(exc))

    if kind == "molecule":
        label = payload.get("label")
        if label in existing["molecule"]:
            errors.append(f"molecule label already exists: {label}")
        parent = payload.get("isotopologue_of")
        if parent and parent not in existing["molecule"] and parent not in staged["molecule"]:
            errors.append(f"unknown isotopologue parent: {parent}")

    if kind == "source":
        nick = payload.get("nick")
        if nick in existing["source"]:
            errors.append(f"source nick already exists: {nick}")

    if kind == "telescope":
        nick = payload.get("nick")
        if nick in existing["telescope"]:
            errors.append(f"telescope nick already exists: {nick}")

    if kind == "detection":
        detection_id = payload.get("id")
        if detection_id in existing["detection"]:
            errors.append(f"detection id already exists: {detection_id}")
        molecule = payload.get("molecule")
        if molecule not in existing["molecule"] and molecule not in staged["molecule"]:
            errors.append(f"unknown molecule: {molecule}")
        for source in payload.get("sources") or []:
            if source not in existing["source"] and source not in staged["source"]:
                errors.append(f"unknown source: {source}")
        for telescope in payload.get("telescopes") or []:
            if telescope not in existing["telescope"] and telescope not in staged["telescope"]:
                errors.append(f"unknown telescope: {telescope}")
        for field_name in DETECTION_RELATION_FIELDS:
            for detection_id in payload.get(field_name) or []:
                if (
                    detection_id not in existing["detection"]
                    and detection_id not in staged["detection"]
                ):
                    errors.append(
                        f"unknown detection relationship target in "
                        f"`{field_name}`: {detection_id}"
                    )

    if kind in REF_ROLE_FIELDS:
        refs = payload.get("refs") or {}
        for role in refs:
            if role not in REF_ROLE_FIELDS[kind]:
                errors.append(f"unknown ref role `{role}`")
        for key in ref_keys_from(payload):
            if key not in ref_keys:
                errors.append(f"unknown reference key: {key}")

    return errors


def identity(kind: str, payload: dict) -> str:
    if kind == "molecule":
        return payload.get("label")
    if kind in {"source", "telescope"}:
        return payload.get("nick")
    if kind == "detection":
        return payload.get("id")
    return "unknown"


def make_report(
    name: str,
    rows: list[dict],
    errors: list[str],
    output_paths: dict[str, Path],
    applied: bool,
) -> str:
    counts = Counter(row["kind"] for row in rows)
    status_counts = Counter(
        row["record"].get("status")
        for row in rows
        if row["kind"] == "detection"
    )
    meta_rows = [row for row in rows if row["meta"]]

    lines = []
    lines.append(f"# Staged Records Preview: {name}")
    lines.append("")
    lines.append("Generated by `scripts/stage_records.py`.")
    lines.append("")
    lines.append("## Summary")
    lines.append("")
    lines.append(f"- Applied to production JSON: `{str(applied).lower()}`")
    lines.append(f"- Staged molecules: {counts['molecule']}")
    lines.append(f"- Staged detections: {counts['detection']}")
    lines.append(f"- Staged sources: {counts['source']}")
    lines.append(f"- Staged telescopes: {counts['telescope']}")
    if status_counts:
        lines.append(
            "- Detection statuses: "
            + ", ".join(f"`{status}`={count}" for status, count in sorted(status_counts.items()))
        )
    lines.append(f"- Validation errors: {len(errors)}")
    lines.append("")
    lines.append("Preview artifacts:")
    for kind, path in output_paths.items():
        lines.append(f"- `{path.relative_to(ROOT)}`")
    lines.append("")

    lines.append("## Validation Errors")
    lines.append("")
    if not errors:
        lines.append("No validation errors.")
    else:
        for error in errors:
            lines.append(f"- {error}")
    lines.append("")

    lines.append("## Staged Records")
    lines.append("")
    for row in rows:
        lines.append(f"### `{row['kind']}` `{row['id']}`")
        lines.append("")
        lines.append(f"- file: `{row['file']}`")
        lines.append(f"- record index: `{row['index']}`")
        if row["kind"] == "detection":
            record = row["record"]
            lines.append(f"- molecule: `{record.get('molecule')}`")
            lines.append(f"- type: `{record.get('type')}`")
            lines.append(f"- status: `{record.get('status')}`")
            lines.append(f"- year: `{record.get('year')}`")
            lines.append(f"- sources: `{', '.join(record.get('sources') or [])}`")
            lines.append(f"- telescopes: `{', '.join(record.get('telescopes') or [])}`")
            lines.append(f"- refs: `{record.get('refs')}`")
            relationships = {
                field: record.get(field)
                for field in DETECTION_RELATION_FIELDS
                if record.get(field)
            }
            if relationships:
                lines.append(f"- relationships: `{relationships}`")
        elif row["kind"] == "molecule":
            record = row["record"]
            lines.append(f"- formula: `{record.get('formula')}`")
            lines.append(f"- name: `{record.get('name')}`")
            lines.append(f"- refs: `{record.get('refs')}`")
        elif row["kind"] == "source":
            record = row["record"]
            lines.append(f"- name: `{record.get('name')}`")
            lines.append(f"- type: `{record.get('type')}`")
        elif row["kind"] == "telescope":
            record = row["record"]
            lines.append(f"- name: `{record.get('name')}`")
            lines.append(f"- type: `{record.get('type')}`")
            lines.append(f"- wavelength: `{', '.join(record.get('wavelength') or [])}`")
        if row["errors"]:
            lines.append(f"- errors: `{'; '.join(row['errors'])}`")
        if row["meta"]:
            lines.append(f"- staging-only notes: `{row['meta']}`")
        lines.append("")

    lines.append("## Staging-Only Notes")
    lines.append("")
    if not meta_rows:
        lines.append("No staging-only notes found.")
    else:
        for row in meta_rows:
            lines.append(f"- `{row['id']}`: `{row['meta']}`")
    lines.append("")

    return "\n".join(lines).rstrip() + "\n"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--staging",
        nargs="+",
        type=Path,
        required=True,
        help="YAML staging file(s) to process.",
    )
    parser.add_argument(
        "--name",
        default=None,
        help="Output name stem. Defaults to the first staging file stem.",
    )
    parser.add_argument(
        "--apply",
        action="store_true",
        help="Write validated staged records into production JSON files.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    name = args.name or args.staging[0].stem

    production = {
        kind: load_json(filename)
        for kind, filename in KIND_TO_FILE.items()
    }
    preview = {kind: json.loads(json.dumps(records)) for kind, records in production.items()}

    existing = {
        "molecule": {record["label"] for record in production["molecule"]},
        "source": {record["nick"] for record in production["source"]},
        "telescope": {record["nick"] for record in production["telescope"]},
    }
    existing["detection"] = {record["id"] for record in production["detection"]}

    staged = defaultdict(set)
    ref_keys = load_reference_keys()

    raw_records = load_yaml_records(args.staging)
    rows = []
    errors = []
    run_date = date.today().isoformat()

    for raw in raw_records:
        kind = raw.get("kind")
        location = f"{raw.get('_staging_file')} record {raw.get('_staging_index')}"
        if kind not in KIND_TO_FILE:
            message = f"{location}: unknown kind `{kind}`"
            errors.append(message)
            rows.append(
                {
                    "kind": kind or "unknown",
                    "id": "unknown",
                    "file": raw.get("_staging_file"),
                    "index": raw.get("_staging_index"),
                    "record": {},
                    "meta": staging_meta(raw),
                    "errors": [message],
                }
            )
            continue

        record = merge_defaults(kind, raw)
        record["history"] = normalize_history(kind, record.get("history"), run_date)
        row_id = identity(kind, record)
        if kind in {"molecule", "source", "telescope", "detection"} and row_id:
            staged[kind].add(row_id)

        rows.append(
            {
                "kind": kind,
                "id": row_id,
                "file": raw.get("_staging_file"),
                "index": raw.get("_staging_index"),
                "record": record,
                "meta": staging_meta(raw),
                "errors": [],
            }
        )

    staged_id_counts = defaultdict(Counter)
    for row in rows:
        if row["kind"] in {"molecule", "source", "telescope", "detection"}:
            staged_id_counts[row["kind"]][row["id"]] += 1

    for row in rows:
        kind = row["kind"]
        if kind not in KIND_TO_FILE:
            continue
        location = f"{row['file']} record {row['index']}"
        row_errors = validate_record(kind, row["record"], existing, staged, ref_keys)
        if (
            kind in {"molecule", "source", "telescope", "detection"}
            and staged_id_counts[kind][row["id"]] > 1
        ):
            row_errors.append(f"duplicate staged {kind} identity: {row['id']}")
        row["errors"] = row_errors
        if row_errors:
            errors.extend(f"{location} `{row['id']}`: {error}" for error in row_errors)
        else:
            preview[kind].append(row["record"])

    output_paths = {
        "molecules": DATA / f"molecules.{name}.preview.json",
        "detections": DATA / f"detections.{name}.preview.json",
        "sources": DATA / f"sources.{name}.preview.json",
        "telescopes": DATA / f"telescopes.{name}.preview.json",
    }
    detail_path = DATA / f"{name}_stage.preview.json"
    report_path = DATA / f"{name}_stage_report.md"

    write_json(output_paths["molecules"], preview["molecule"])
    write_json(output_paths["detections"], preview["detection"])
    write_json(output_paths["sources"], preview["source"])
    write_json(output_paths["telescopes"], preview["telescope"])
    write_json(detail_path, rows)

    report_paths = dict(output_paths)
    report_paths["details"] = detail_path
    report_path.write_text(make_report(name, rows, errors, report_paths, applied=args.apply))

    if errors:
        print(f"Wrote {report_path.relative_to(ROOT)} with {len(errors)} validation errors.")
        raise SystemExit(1)

    if args.apply:
        staged_kinds = {
            row["kind"]
            for row in rows
            if row["kind"] in KIND_TO_FILE and not row["errors"]
        }
        for kind in staged_kinds:
            write_json(DATA / KIND_TO_FILE[kind], preview[kind])

    print(f"Wrote {report_path.relative_to(ROOT)}")
    for path in [*output_paths.values(), detail_path]:
        print(f"Wrote {path.relative_to(ROOT)}")
    if args.apply:
        print("Applied staged records to production JSON.")


if __name__ == "__main__":
    main()
