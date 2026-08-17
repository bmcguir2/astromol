from __future__ import annotations

import importlib.util
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]
TEMPLATE_DIR = ROOT / "curation" / "templates"
STAGE_RECORDS_PATH = ROOT / "scripts" / "stage_records.py"


def load_stage_records_module():
    spec = importlib.util.spec_from_file_location(
        "stage_records",
        STAGE_RECORDS_PATH,
    )
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_template_record(kind: str) -> dict:
    records = yaml.safe_load((TEMPLATE_DIR / f"{kind}.yaml").read_text())
    assert isinstance(records, list)
    assert len(records) == 1
    record = records[0]
    assert isinstance(record, dict)
    assert record["kind"] == kind
    return record


def test_templates_match_staging_fields():
    stage_records = load_stage_records_module()

    for kind, expected_fields in stage_records.FIELDS.items():
        record = load_template_record(kind)
        assert record["operation"] == "add"
        template_fields = {
            key
            for key in record
            if key not in stage_records.STAGING_CONTROL_FIELDS
            and not key.startswith("_")
        }

        assert template_fields == set(expected_fields)


def test_template_extra_fields_are_staging_only():
    stage_records = load_stage_records_module()

    for kind in stage_records.FIELDS:
        record = load_template_record(kind)
        extra_fields = [
            key
            for key in record
            if key not in stage_records.STAGING_CONTROL_FIELDS
            and key not in stage_records.FIELDS[kind]
        ]

        assert all(key.startswith("_") for key in extra_fields)


def test_required_and_default_fields_are_known_staging_fields():
    stage_records = load_stage_records_module()

    for kind, required_fields in stage_records.REQUIRED.items():
        assert set(required_fields).issubset(stage_records.FIELDS[kind])

    for kind, defaults in stage_records.DEFAULTS.items():
        assert set(defaults).issubset(stage_records.FIELDS[kind])
