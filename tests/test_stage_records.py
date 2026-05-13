from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
STAGE_RECORDS_PATH = ROOT / "scripts" / "stage_records.py"


def load_stage_records_module():
    spec = importlib.util.spec_from_file_location(
        "stage_records_for_tests",
        STAGE_RECORDS_PATH,
    )
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_json(path: Path, payload):
    path.write_text(json.dumps(payload, indent=2) + "\n")


def make_temp_data_dir(tmp_path: Path) -> Path:
    data = tmp_path / "astromol" / "data"
    data.mkdir(parents=True)

    write_json(data / "molecules.json", [])
    write_json(data / "detections.json", [])
    write_json(data / "sources.json", [])
    write_json(data / "telescopes.json", [])
    (data / "references.bib").write_text(
        "@article{Example:2026:1,\n"
        "  author = {Example, A.},\n"
        "  year = {2026},\n"
        "  journal = {ApJ}\n"
        "}\n"
    )
    return data


def run_stage_records(monkeypatch, tmp_path: Path, yaml_text: str, *args: str):
    stage_records = load_stage_records_module()
    data = make_temp_data_dir(tmp_path)
    staging_file = tmp_path / "example.yaml"
    staging_file.write_text(yaml_text)

    monkeypatch.setattr(stage_records, "ROOT", tmp_path)
    monkeypatch.setattr(stage_records, "DATA", data)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "stage_records.py",
            "--staging",
            str(staging_file),
            "--name",
            "example",
            *args,
        ],
    )

    stage_records.main()
    return data


def valid_staging_yaml() -> str:
    return """
- kind: source
  name: Example Source
  nick: ExampleSource
  type: Dark Cloud
  ra: "12:34:56"
  dec: "-01:23:45"

- kind: telescope
  name: Example Telescope
  nick: ExampleTelescope
  shortname: Example 12-m
  type: Single Dish
  wavelength: [mm]

- kind: molecule
  label: mol:EXAMPLE
  name: example molecule
  formula: CH
  table_formula:
  refs:
    lab: [Example:2026:1]
  _curation_note: keep this out of production JSON

- kind: detection
  id: det:EXAMPLE:ism-csm:2026
  molecule: mol:EXAMPLE
  sources: [ExampleSource]
  telescopes: [ExampleTelescope]
  wavelengths: [mm]
  year: 2026
  type: ISM/CSM
  refs:
    observation: [Example:2026:1]
"""


def test_stage_records_writes_preview_without_modifying_production(
    monkeypatch,
    tmp_path,
):
    data = run_stage_records(monkeypatch, tmp_path, valid_staging_yaml())

    assert json.loads((data / "molecules.json").read_text()) == []
    assert json.loads((data / "detections.json").read_text()) == []
    assert json.loads((data / "sources.json").read_text()) == []
    assert json.loads((data / "telescopes.json").read_text()) == []

    molecule_preview = json.loads(
        (data / "molecules.example.preview.json").read_text()
    )
    detection_preview = json.loads(
        (data / "detections.example.preview.json").read_text()
    )
    source_preview = json.loads((data / "sources.example.preview.json").read_text())
    telescope_preview = json.loads(
        (data / "telescopes.example.preview.json").read_text()
    )

    assert molecule_preview[0]["label"] == "mol:EXAMPLE"
    assert molecule_preview[0]["table_formula"] == "CH"
    assert molecule_preview[0]["history"]["accepted"]["census"] == "2026"
    assert "_curation_note" not in molecule_preview[0]

    assert detection_preview[0]["id"] == "det:EXAMPLE:ism-csm:2026"
    assert detection_preview[0]["history"]["accepted"]["context"] == "confirmed_ism_csm"
    assert source_preview[0]["nick"] == "ExampleSource"
    assert "accepted" not in source_preview[0]["history"]
    assert telescope_preview[0]["nick"] == "ExampleTelescope"
    assert "accepted" not in telescope_preview[0]["history"]

    detail = json.loads((data / "example_stage.preview.json").read_text())
    molecule_row = next(row for row in detail if row["id"] == "mol:EXAMPLE")
    assert molecule_row["meta"] == {
        "_curation_note": "keep this out of production JSON"
    }

    report = (data / "example_stage_report.md").read_text()
    assert "Validation errors: 0" in report
    assert "Applied to production JSON: `false`" in report


def test_stage_records_apply_writes_valid_records_to_production(
    monkeypatch,
    tmp_path,
):
    data = run_stage_records(monkeypatch, tmp_path, valid_staging_yaml(), "--apply")

    assert json.loads((data / "molecules.json").read_text())[0]["label"] == "mol:EXAMPLE"
    assert json.loads((data / "detections.json").read_text())[0]["id"] == (
        "det:EXAMPLE:ism-csm:2026"
    )
    assert json.loads((data / "sources.json").read_text())[0]["nick"] == "ExampleSource"
    assert json.loads((data / "telescopes.json").read_text())[0]["nick"] == (
        "ExampleTelescope"
    )

    report = (data / "example_stage_report.md").read_text()
    assert "Applied to production JSON: `true`" in report


def test_stage_records_rejects_unknown_references(monkeypatch, tmp_path):
    yaml_text = """
- kind: molecule
  label: mol:BADREF
  name: bad reference molecule
  formula: CH
  refs:
    lab: [Missing:2026:1]
"""
    stage_records = load_stage_records_module()
    data = make_temp_data_dir(tmp_path)
    staging_file = tmp_path / "badref.yaml"
    staging_file.write_text(yaml_text)

    monkeypatch.setattr(stage_records, "ROOT", tmp_path)
    monkeypatch.setattr(stage_records, "DATA", data)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "stage_records.py",
            "--staging",
            str(staging_file),
            "--name",
            "badref",
        ],
    )

    try:
        stage_records.main()
    except SystemExit as exc:
        assert exc.code == 1
    else:  # pragma: no cover - explicit failure branch
        raise AssertionError("Expected staging validation to fail")

    assert json.loads((data / "molecules.json").read_text()) == []
    report = (data / "badref_stage_report.md").read_text()
    assert "unknown reference key: Missing:2026:1" in report

