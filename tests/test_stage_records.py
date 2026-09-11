from __future__ import annotations

from datetime import date
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import pytest
import yaml


ROOT = Path(__file__).resolve().parents[1]
STAGE_RECORDS_PATH = ROOT / "scripts" / "stage_records.py"
CLEANUP_STAGE_PATH = ROOT / "scripts" / "cleanup_stage.py"
TEST_RUN_DATE = date(2026, 8, 17)


class FixedDate(date):
    @classmethod
    def today(cls):
        return cls.fromordinal(TEST_RUN_DATE.toordinal())


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


def run_stage_records(
    monkeypatch,
    tmp_path: Path,
    yaml_text: str,
    *args: str,
    data: Path | None = None,
):
    stage_records = load_stage_records_module()
    data = data or make_temp_data_dir(tmp_path)
    staging_file = tmp_path / "example.yaml"
    staging_file.write_text(yaml_text)

    monkeypatch.setattr(stage_records, "ROOT", tmp_path)
    monkeypatch.setattr(stage_records, "DATA", data)
    monkeypatch.setattr(stage_records, "date", FixedDate)
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
    assert set(manifest["production_hashes"]) == set(manifest["production_files"])

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


def test_prepare_update_writes_full_locked_template(monkeypatch, tmp_path):
    data = run_stage_records(
        monkeypatch,
        tmp_path,
        valid_staging_yaml(),
        "--apply",
    )
    stage_records = load_stage_records_module()
    monkeypatch.setattr(stage_records, "ROOT", tmp_path)
    monkeypatch.setattr(stage_records, "DATA", data)
    output = tmp_path / "molecule_update.yaml"

    stage_records.prepare_update_template("molecule", "mol:EXAMPLE", output)

    [record] = yaml.safe_load(output.read_text())
    production_record = json.loads((data / "molecules.json").read_text())[0]
    assert record["kind"] == "molecule"
    assert record["operation"] == "update"
    assert record["_base_digest"] == stage_records.record_digest(production_record)
    assert record["_event_kind"] == "updated"
    assert record["_update_summary"] == ""
    assert set(stage_records.FIELDS["molecule"]).issubset(record)


def test_stage_records_updates_existing_full_record(monkeypatch, tmp_path):
    data = run_stage_records(
        monkeypatch,
        tmp_path,
        valid_staging_yaml(),
        "--apply",
    )
    stage_records = load_stage_records_module()
    molecule = json.loads((data / "molecules.json").read_text())[0]
    update = {
        "kind": "molecule",
        "operation": "update",
        "_base_digest": stage_records.record_digest(molecule),
        "_event_kind": "updated",
        "_update_summary": "Add curator note.",
        **molecule,
    }
    update["note"] = "Updated through the staging workflow."
    yaml_text = yaml.safe_dump([update], sort_keys=False)

    run_stage_records(monkeypatch, tmp_path, yaml_text, data=data)

    assert json.loads((data / "molecules.json").read_text())[0]["note"] is None
    preview = json.loads((data / "molecules.example.preview.json").read_text())[0]
    assert preview["note"] == "Updated through the staging workflow."
    assert preview["history"]["events"][-1] == {
        "kind": "updated",
        "summary": "Add curator note.",
        "date": "2026-08-17",
        "fields": ["note"],
    }
    report = (data / "example_stage_report.md").read_text()
    assert "operation: `update`" in report
    assert "changed fields: `note`" in report

    run_stage_records(monkeypatch, tmp_path, yaml_text, "--apply", data=data)
    applied = json.loads((data / "molecules.json").read_text())[0]
    assert applied["note"] == "Updated through the staging workflow."
    manifest = json.loads((data / "example_stage_manifest.json").read_text())
    assert "astromol/data/molecules.json" in manifest["production_hashes"]


def test_stage_records_rejects_stale_update_digest(monkeypatch, tmp_path):
    data = run_stage_records(
        monkeypatch,
        tmp_path,
        valid_staging_yaml(),
        "--apply",
    )
    molecule = json.loads((data / "molecules.json").read_text())[0]
    update = {
        "kind": "molecule",
        "operation": "update",
        "_base_digest": "stale",
        "_event_kind": "updated",
        "_update_summary": "Stale update.",
        **molecule,
    }
    update["note"] = "This must not apply."

    with pytest.raises(SystemExit):
        run_stage_records(
            monkeypatch,
            tmp_path,
            yaml.safe_dump([update], sort_keys=False),
            "--apply",
            data=data,
        )

    assert json.loads((data / "molecules.json").read_text())[0]["note"] is None
    report = (data / "example_stage_report.md").read_text()
    assert "stale `_base_digest`" in report


def test_stage_records_promotes_existing_molecule_explicitly(monkeypatch, tmp_path):
    initial_yaml = """
- kind: molecule
  label: mol:EXAMPLE
  name: example molecule
  formula: CH
  history:
    introduced:
      context: tentative
    accepted: null
"""
    data = run_stage_records(monkeypatch, tmp_path, initial_yaml, "--apply")
    stage_records = load_stage_records_module()
    molecule = json.loads((data / "molecules.json").read_text())[0]
    update = {
        "kind": "molecule",
        "operation": "update",
        "_base_digest": stage_records.record_digest(molecule),
        "_event_kind": "updated",
        "_update_summary": "Promote molecule after secure ice detection.",
        **molecule,
    }
    update["history"] = json.loads(json.dumps(molecule["history"]))
    update["history"]["accepted"] = {
        "date": "2026-08-17",
        "census": "2026",
        "context": "confirmed_ice",
    }

    run_stage_records(
        monkeypatch,
        tmp_path,
        yaml.safe_dump([update], sort_keys=False),
        "--apply",
        data=data,
    )

    promoted = json.loads((data / "molecules.json").read_text())[0]
    assert promoted["history"]["accepted"]["context"] == "confirmed_ice"
    assert promoted["history"]["events"][-1]["fields"] == [
        "history.accepted"
    ]


def test_stage_records_derives_reciprocal_confirmation(monkeypatch, tmp_path):
    initial_yaml = """
- kind: molecule
  label: mol:EXAMPLE
  name: example molecule
  formula: CH

- kind: detection
  id: det:EXAMPLE:ice:2005
  molecule: mol:EXAMPLE
  year: 2005
  type: ice
  status: tentative
  refs:
    observation: [Example:2026:1]
"""
    data = run_stage_records(monkeypatch, tmp_path, initial_yaml, "--apply")
    confirmation_yaml = """
- kind: detection
  id: det:EXAMPLE:ice:2024
  molecule: mol:EXAMPLE
  year: 2024
  type: ice
  status: secure
  first: true
  refs:
    observation: [Example:2026:1]
  confirms: [det:EXAMPLE:ice:2005]
"""

    run_stage_records(monkeypatch, tmp_path, confirmation_yaml, data=data)

    preview = {
        record["id"]: record
        for record in json.loads(
            (data / "detections.example.preview.json").read_text()
        )
    }
    assert preview["det:EXAMPLE:ice:2024"]["confirms"] == [
        "det:EXAMPLE:ice:2005"
    ]
    assert preview["det:EXAMPLE:ice:2005"]["confirmed_by"] == [
        "det:EXAMPLE:ice:2024"
    ]
    assert preview["det:EXAMPLE:ice:2005"]["history"]["events"][-1][
        "fields"
    ] == ["confirmed_by"]
    report = (data / "example_stage_report.md").read_text()
    assert "Derived reciprocal updates: 1" in report
    assert "`det:EXAMPLE:ice:2005.confirmed_by`: `add`" in report

    run_stage_records(
        monkeypatch,
        tmp_path,
        confirmation_yaml,
        "--apply",
        data=data,
    )
    applied = {
        record["id"]: record
        for record in json.loads((data / "detections.json").read_text())
    }
    assert applied["det:EXAMPLE:ice:2005"]["confirmed_by"] == [
        "det:EXAMPLE:ice:2024"
    ]


@pytest.mark.parametrize(
    ("field_name", "reciprocal_name"),
    [
        ("confirms", "confirmed_by"),
        ("confirmed_by", "confirms"),
        ("disputes", "disputed_by"),
        ("disputed_by", "disputes"),
        ("supersedes", "superseded_by"),
        ("superseded_by", "supersedes"),
    ],
)
def test_all_detection_relationships_derive_reciprocals(
    field_name,
    reciprocal_name,
):
    stage_records = load_stage_records_module()
    target = {
        "id": "det:TARGET:ice:2000",
        reciprocal_name: [],
        "history": {"events": []},
    }
    source = {
        "id": "det:SOURCE:ice:2020",
        field_name: [target["id"]],
        "history": {"events": []},
    }
    production = {
        "molecule": [],
        "source": [],
        "telescope": [],
        "detection": [target],
    }
    preview = {
        **production,
        "detection": [json.loads(json.dumps(target)), source],
    }
    rows = [
        {
            "kind": "detection",
            "errors": [],
            "before": None,
            "record": source,
        }
    ]

    derived = stage_records.apply_derived_reciprocals(
        preview,
        rows,
        production,
        "2026-08-17",
    )

    assert preview["detection"][0][reciprocal_name] == [source["id"]]
    assert derived == [
        {
            "record": target["id"],
            "field": reciprocal_name,
            "action": "add",
            "related_record": source["id"],
        }
    ]

    before_source = json.loads(json.dumps(source))
    after_source = json.loads(json.dumps(source))
    after_source[field_name] = []
    remove_rows = [
        {
            "kind": "detection",
            "errors": [],
            "before": before_source,
            "record": after_source,
        }
    ]
    removed = stage_records.apply_derived_reciprocals(
        preview,
        remove_rows,
        production,
        "2026-08-17",
    )
    assert preview["detection"][0][reciprocal_name] == []
    assert removed == [
        {
            "record": target["id"],
            "field": reciprocal_name,
            "action": "remove",
            "related_record": source["id"],
        }
    ]


def test_semantic_preview_errors_block_apply(monkeypatch, tmp_path):
    yaml_text = """
- kind: molecule
  label: mol:EXAMPLE
  name: example molecule
  formula: CH

- kind: detection
  id: det:EXAMPLE:ice:2025
  molecule: mol:EXAMPLE
  year: 2025
  type: ice
  first: true
  refs:
    observation: [Example:2026:1]

- kind: detection
  id: det:EXAMPLE:ice:2026
  molecule: mol:EXAMPLE
  year: 2026
  type: ice
  first: true
  refs:
    observation: [Example:2026:1]
"""
    data = make_temp_data_dir(tmp_path)

    with pytest.raises(SystemExit):
        run_stage_records(
            monkeypatch,
            tmp_path,
            yaml_text,
            "--apply",
            data=data,
        )

    assert json.loads((data / "molecules.json").read_text()) == []
    assert json.loads((data / "detections.json").read_text()) == []
    assert not (data / "example_stage_manifest.json").exists()
    report = (data / "example_stage_report.md").read_text()
    assert "semantic validation `first-flag-extra`" in report


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


def test_cleanup_stage_verification_failure_preserves_artifacts(
    monkeypatch,
    tmp_path,
):
    data = run_stage_records(monkeypatch, tmp_path, valid_staging_yaml(), "--apply")
    cleanup_stage = load_cleanup_stage_module()
    manifest_path = data / "example_stage_manifest.json"
    staging_file = tmp_path / "example.yaml"

    monkeypatch.setattr(cleanup_stage, "ROOT", tmp_path)
    monkeypatch.setattr(cleanup_stage, "DATA", data)
    monkeypatch.setattr(
        cleanup_stage,
        "run_curation_verification",
        lambda: (_ for _ in ()).throw(RuntimeError("verification failed")),
    )
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "cleanup_stage.py",
            "--name",
            "example",
            "--commit-message",
            "Do not commit",
        ],
    )

    with pytest.raises(RuntimeError, match="verification failed"):
        cleanup_stage.main()

    assert manifest_path.exists()
    assert staging_file.exists()
    assert (data / "example_stage_report.md").exists()


def test_cleanup_stage_can_commit_deleted_and_curated_files(monkeypatch, tmp_path):
    data = run_stage_records(monkeypatch, tmp_path, valid_staging_yaml(), "--apply")
    cleanup_stage = load_cleanup_stage_module()
    verification_runs = []
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
        cleanup_stage,
        "run_curation_verification",
        lambda: verification_runs.append(True),
    )
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
    assert verification_runs == [True]

    status = subprocess.run(
        ["git", "status", "--short"],
        cwd=tmp_path,
        check=True,
        capture_output=True,
        text=True,
    )
    assert status.stdout.strip() == ""


def test_cleanup_stage_rejects_modified_applied_file(monkeypatch, tmp_path):
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
            "--skip-verification",
        ],
    )

    with pytest.raises(SystemExit, match="content changed after apply"):
        cleanup_stage.main()

    assert (data / "example_stage_manifest.json").exists()
    assert (tmp_path / "example.yaml").exists()
