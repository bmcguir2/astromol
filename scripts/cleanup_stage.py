"""Clean up staging YAML and preview artifacts from a stage_records run.

This script reads a manifest written by scripts/stage_records.py, deletes the
listed staging and preview files, then deletes the manifest itself. It can also
optionally stage the related git paths, create a commit, and push the current
branch.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "astromol" / "data"
DEFAULT_COMMIT_INCLUDE_PATHS = [
    Path("astromol/data/references.bib"),
    Path("tests/baselines/production_data.json"),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    location = parser.add_mutually_exclusive_group(required=True)
    location.add_argument(
        "--manifest",
        type=Path,
        help="Manifest JSON written by scripts/stage_records.py.",
    )
    location.add_argument(
        "--name",
        help="Output name stem used by scripts/stage_records.py.",
    )
    parser.add_argument(
        "--include",
        action="append",
        type=Path,
        default=[],
        help="Additional repo-relative paths to stage if committing.",
    )
    parser.add_argument(
        "--commit-message",
        default=None,
        help="Commit staged cleanup and curation files with this message.",
    )
    parser.add_argument(
        "--push",
        action="store_true",
        help="Push after a successful commit. Requires --commit-message.",
    )
    args = parser.parse_args()
    if args.push and not args.commit_message:
        parser.error("--push requires --commit-message")
    return args


def manifest_path_from_args(args: argparse.Namespace) -> Path:
    if args.manifest is not None:
        return args.manifest
    return DATA / f"{args.name}_stage_manifest.json"


def load_manifest(path: Path) -> dict:
    return json.loads(path.read_text())


def resolve_repo_path(path_text: str) -> Path:
    path = (ROOT / path_text).resolve()
    root = ROOT.resolve()
    try:
        path.relative_to(root)
    except ValueError as exc:
        raise ValueError(f"manifest path escapes repository root: {path_text}") from exc
    return path


def git_status_entries(path: Path) -> str:
    result = subprocess.run(
        ["git", "status", "--short", "--", str(path.relative_to(ROOT))],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip()


def git_add_paths(paths: list[Path]) -> list[Path]:
    staged = []
    seen = set()
    for path in paths:
        resolved = path.resolve()
        if resolved in seen:
            continue
        seen.add(resolved)
        if not git_status_entries(path):
            continue
        subprocess.run(
            ["git", "add", "--", str(path.relative_to(ROOT))],
            cwd=ROOT,
            check=True,
        )
        staged.append(path)
    return staged


def main() -> None:
    args = parse_args()
    manifest_path = manifest_path_from_args(args).resolve()
    manifest = load_manifest(manifest_path)

    delete_paths = [
        resolve_repo_path(path_text)
        for path_text in (
            list(manifest.get("staging_files", []))
            + list(manifest.get("preview_artifacts", []))
        )
    ]

    deleted = []
    for path in delete_paths:
        if path.exists():
            path.unlink()
            deleted.append(path)

    if manifest_path.exists():
        manifest_path.unlink()
        deleted.append(manifest_path)

    print(f"Deleted {len(deleted)} staging file(s).")
    for path in deleted:
        print(f"Deleted {path.relative_to(ROOT)}")

    if not args.commit_message:
        return

    stage_candidates = [
        *deleted,
        *[
            resolve_repo_path(path_text)
            for path_text in manifest.get("production_files", [])
        ],
        *[(ROOT / path).resolve() for path in DEFAULT_COMMIT_INCLUDE_PATHS],
        *[resolve_repo_path(str(path)) for path in args.include],
    ]
    staged = git_add_paths(stage_candidates)
    if not staged:
        raise SystemExit("No git changes matched the cleanup manifest.")

    subprocess.run(
        ["git", "commit", "-m", args.commit_message],
        cwd=ROOT,
        check=True,
    )
    print(f"Created commit: {args.commit_message}")

    if args.push:
        subprocess.run(
            ["git", "push"],
            cwd=ROOT,
            check=True,
        )
        print("Pushed current branch.")


if __name__ == "__main__":
    main()
