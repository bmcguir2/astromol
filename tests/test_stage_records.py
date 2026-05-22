from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import pytest


ROOT = Path(__file__).resolve().parents[1]
STAGE_RECORDS_PATH = ROOT / "scripts" / "stage_records.py"
CLEANUP_STAGE_PATH = ROOT / "scripts" / "cleanup_stage.py"


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


def load_cleanup_stage_module():
    spec = importlib.util.spec_from_file_location(
        "cleanup_stage_for_tests",
        CLEANUP_STAGE_PATH,
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
    assert "## Generated Count Updates" in report
    assert "`counts.molecules`: `0` -> `1`" in report
    assert "`regression_counts.census_view_2026.ism_molecules`: `0` -> `1`" in report

    manifest = json.loads((data / "example_stage_manifest.json").read_text())
    assert manifest["name"] == "example"
    assert manifest["applied"] is False
    assert manifest["staging_files"] == ["example.yaml"]
    assert "astromol/data/example_stage_report.md" in manifest["preview_artifacts"]
    assert manifest["production_files"] == []


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

    manifest = json.loads((data / "example_stage_manifest.json").read_text())
    assert manifest["applied"] is True
    assert sorted(manifest["production_files"]) == [
        "astromol/data/detections.json",
        "astromol/data/molecules.json",
        "astromol/data/sources.json",
        "astromol/data/telescopes.json",
        "tests/baselines/production_data.json",
    ]

    baseline = json.loads(
        (tmp_path / "tests" / "baselines" / "production_data.json").read_text()
    )
    assert baseline["counts"]["molecules"] == 1
    assert baseline["regression_counts"]["census_view_2026"]["ism_molecules"] == 1


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
    assert not (data / "badref_stage_manifest.json").exists()


def init_git_repo(path: Path) -> None:
    subprocess.run(["git", "init"], cwd=path, check=True, capture_output=True, text=True)
    subprocess.run(
        ["git", "config", "user.name", "Test User"],
        cwd=path,
        check=True,
        capture_output=True,
        text=True,
    )
    subprocess.run(
        ["git", "config", "user.email", "test@example.com"],
        cwd=path,
        check=True,
        capture_output=True,
        text=True,
    )


def test_cleanup_stage_parses_github_remote_urls():
    cleanup_stage = load_cleanup_stage_module()

    assert (
        cleanup_stage.github_repo_from_remote_url(
            "https://github.com/bmcguir2/astromol.git"
        )
        == "bmcguir2/astromol"
    )
    assert (
        cleanup_stage.github_repo_from_remote_url(
            "git@github.com:bmcguir2/astromol.git"
        )
        == "bmcguir2/astromol"
    )
    assert (
        cleanup_stage.github_repo_from_remote_url(
            "ssh://git@github.com/bmcguir2/astromol.git"
        )
        == "bmcguir2/astromol"
    )


def test_cleanup_stage_close_issue_requires_push(monkeypatch):
    cleanup_stage = load_cleanup_stage_module()
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "cleanup_stage.py",
            "--name",
            "example",
            "--close-issue",
            "123",
        ],
    )

    with pytest.raises(SystemExit):
        cleanup_stage.parse_args()


def test_cleanup_stage_closes_github_issue_with_commit_link(monkeypatch, tmp_path):
    cleanup_stage = load_cleanup_stage_module()
    commands = []

    def fake_run(args, **kwargs):
        commands.append(args)
        return subprocess.CompletedProcess(args, 0, stdout="", stderr="")

    monkeypatch.setattr(cleanup_stage, "ROOT", tmp_path)
    monkeypatch.setattr(
        cleanup_stage,
        "git_output",
        lambda args: "https://github.com/bmcguir2/astromol.git",
    )
    monkeypatch.setattr(cleanup_stage.shutil, "which", lambda name: "/usr/bin/gh")
    monkeypatch.setattr(cleanup_stage.subprocess, "run", fake_run)

    cleanup_stage.close_github_issues(
        [123],
        "2026_thioacetaldehyde_tmc1",
        "abcdef1234567890",
    )

    assert commands == [
        [
            "gh",
            "issue",
            "close",
            "123",
            "--repo",
            "bmcguir2/astromol",
            "--reason",
            "completed",
            "--comment",
            (
                "Applied the `2026_thioacetaldehyde_tmc1` curation batch in commit "
                "[`abcdef1`](https://github.com/bmcguir2/astromol/commit/"
                "abcdef1234567890).\n\n"
                "The staged records were applied, committed, and pushed with the "
                "refreshed production-data baseline."
            ),
        ]
    ]


def test_cleanup_stage_removes_manifest_listed_files(monkeypatch, tmp_path):
    data = run_stage_records(monkeypatch, tmp_path, valid_staging_yaml())
    cleanup_stage = load_cleanup_stage_module()

    manifest_path = data / "example_stage_manifest.json"
    staging_file = tmp_path / "example.yaml"
    assert manifest_path.exists()
    assert staging_file.exists()
    assert (data / "molecules.example.preview.json").exists()
    assert (data / "example_stage_report.md").exists()

    monkeypatch.setattr(cleanup_stage, "ROOT", tmp_path)
    monkeypatch.setattr(cleanup_stage, "DATA", data)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "cleanup_stage.py",
            "--name",
            "example",
        ],
    )

    cleanup_stage.main()

    assert not manifest_path.exists()
    assert not staging_file.exists()
    assert not (data / "molecules.example.preview.json").exists()
    assert not (data / "detections.example.preview.json").exists()
    assert not (data / "sources.example.preview.json").exists()
    assert not (data / "telescopes.example.preview.json").exists()
    assert not (data / "example_stage.preview.json").exists()
    assert not (data / "example_stage_report.md").exists()


def test_cleanup_stage_can_commit_deleted_and_curated_files(monkeypatch, tmp_path):
    data = run_stage_records(monkeypatch, tmp_path, valid_staging_yaml(), "--apply")
    cleanup_stage = load_cleanup_stage_module()
    init_git_repo(tmp_path)

    tracked_paths = [
        "astromol/data/molecules.json",
        "astromol/data/detections.json",
        "astromol/data/sources.json",
        "astromol/data/telescopes.json",
        "astromol/data/references.bib",
        "example.yaml",
        "astromol/data/molecules.example.preview.json",
        "astromol/data/detections.example.preview.json",
        "astromol/data/sources.example.preview.json",
        "astromol/data/telescopes.example.preview.json",
        "astromol/data/example_stage.preview.json",
        "astromol/data/example_stage_report.md",
        "astromol/data/example_stage_manifest.json",
    ]
    subprocess.run(
        ["git", "add", "--", *tracked_paths],
        cwd=tmp_path,
        check=True,
        capture_output=True,
        text=True,
    )
    subprocess.run(
        ["git", "commit", "-m", "baseline"],
        cwd=tmp_path,
        check=True,
        capture_output=True,
        text=True,
    )

    (data / "references.bib").write_text(
        (data / "references.bib").read_text() + "\n@article{Second:2026:2,\n  author = {Second, B.}\n}\n"
    )

    monkeypatch.setattr(cleanup_stage, "ROOT", tmp_path)
    monkeypatch.setattr(cleanup_stage, "DATA", data)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "cleanup_stage.py",
            "--name",
            "example",
            "--commit-message",
            "Clean staged example curation files",
        ],
    )

    cleanup_stage.main()

    log = subprocess.run(
        ["git", "log", "-1", "--pretty=%s"],
        cwd=tmp_path,
        check=True,
        capture_output=True,
        text=True,
    )
    assert log.stdout.strip() == "Clean staged example curation files"

    status = subprocess.run(
        ["git", "status", "--short"],
        cwd=tmp_path,
        check=True,
        capture_output=True,
        text=True,
    )
    assert status.stdout.strip() == ""


def test_cleanup_stage_auto_stages_baseline_when_modified(monkeypatch, tmp_path):
    data = run_stage_records(monkeypatch, tmp_path, valid_staging_yaml(), "--apply")
    cleanup_stage = load_cleanup_stage_module()
    init_git_repo(tmp_path)

    baseline_dir = tmp_path / "tests" / "baselines"
    baseline_dir.mkdir(parents=True, exist_ok=True)
    baseline_path = baseline_dir / "production_data.json"
    baseline_path.write_text('{"counts": 1}\n')

    tracked_paths = [
        "astromol/data/molecules.json",
        "astromol/data/detections.json",
        "astromol/data/sources.json",
        "astromol/data/telescopes.json",
        "astromol/data/references.bib",
        "tests/baselines/production_data.json",
        "example.yaml",
        "astromol/data/molecules.example.preview.json",
        "astromol/data/detections.example.preview.json",
        "astromol/data/sources.example.preview.json",
        "astromol/data/telescopes.example.preview.json",
        "astromol/data/example_stage.preview.json",
        "astromol/data/example_stage_report.md",
        "astromol/data/example_stage_manifest.json",
    ]
    subprocess.run(
        ["git", "add", "--", *tracked_paths],
        cwd=tmp_path,
        check=True,
        capture_output=True,
        text=True,
    )
    subprocess.run(
        ["git", "commit", "-m", "baseline"],
        cwd=tmp_path,
        check=True,
        capture_output=True,
        text=True,
    )

    baseline_path.write_text('{"counts": 2}\n')

    monkeypatch.setattr(cleanup_stage, "ROOT", tmp_path)
    monkeypatch.setattr(cleanup_stage, "DATA", data)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "cleanup_stage.py",
            "--name",
            "example",
            "--commit-message",
            "Clean staged example curation files",
        ],
    )

    cleanup_stage.main()

    show = subprocess.run(
        ["git", "show", "--stat", "--oneline", "-1"],
        cwd=tmp_path,
        check=True,
        capture_output=True,
        text=True,
    )
    assert "tests/baselines/production_data.json" in show.stdout
