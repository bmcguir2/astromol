"""Clean up staging YAML and preview artifacts from a stage_records run.

This script reads a manifest written by scripts/stage_records.py, deletes the
listed staging and preview files, then deletes the manifest itself. It can also
optionally stage the related git paths, create a commit, and push the current
branch.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
from urllib.parse import urlparse


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
        "--skip-verification",
        action="store_true",
        help="Skip the default scripts/check_curation.py gate before committing.",
    )
    parser.add_argument(
        "--push",
        action="store_true",
        help="Push after a successful commit. Requires --commit-message.",
    )
    parser.add_argument(
        "--close-issue",
        action="append",
        type=int,
        default=[],
        metavar="NUMBER",
        help=(
            "Close a GitHub issue after a successful push. May be passed more "
            "than once. Requires --push."
        ),
    )
    args = parser.parse_args()
    if args.push and not args.commit_message:
        parser.error("--push requires --commit-message")
    if args.close_issue and not args.push:
        parser.error("--close-issue requires --push")
    return args


def manifest_path_from_args(args: argparse.Namespace) -> Path:
    if args.manifest is not None:
        return args.manifest
    return DATA / f"{args.name}_stage_manifest.json"


def load_manifest(path: Path) -> dict:
    return json.loads(path.read_text())


def file_digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_manifest_hashes(manifest: dict) -> None:
    """Reject production drift after an applied staging batch."""
    mismatches = []
    for path_text, expected_digest in manifest.get("production_hashes", {}).items():
        path = resolve_repo_path(path_text)
        if not path.exists():
            mismatches.append(f"{path_text}: file is missing")
        elif file_digest(path) != expected_digest:
            mismatches.append(f"{path_text}: content changed after apply")
    if mismatches:
        raise SystemExit(
            "Applied production files no longer match the staging manifest:\n- "
            + "\n- ".join(mismatches)
            + "\nRegenerate/reapply the staging batch before cleanup."
        )


def run_curation_verification() -> None:
    command = [sys.executable, "scripts/check_curation.py"]
    subprocess.run(command, cwd=ROOT, check=True)


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


def git_output(args: list[str]) -> str:
    result = subprocess.run(
        ["git", *args],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip()


def stage_name_from_args(args: argparse.Namespace, manifest_path: Path) -> str:
    if args.name:
        return args.name
    stem = manifest_path.name
    suffix = "_stage_manifest.json"
    if stem.endswith(suffix):
        return stem[: -len(suffix)]
    return manifest_path.stem


def github_repo_from_remote_url(remote_url: str) -> str:
    remote_url = remote_url.strip()
    path = ""
    if remote_url.startswith("git@github.com:"):
        path = remote_url.split(":", 1)[1]
    elif remote_url.startswith("ssh://git@github.com/"):
        path = urlparse(remote_url).path.lstrip("/")
    else:
        parsed = urlparse(remote_url)
        if parsed.netloc.lower() == "github.com":
            path = parsed.path.lstrip("/")

    if path.endswith(".git"):
        path = path[:-4]

    parts = [part for part in path.split("/") if part]
    if len(parts) != 2:
        raise ValueError(
            f"cannot determine GitHub repository from origin URL: {remote_url}"
        )
    return "/".join(parts)


def issue_close_comment(stage_name: str, commit_sha: str, commit_url: str) -> str:
    short_sha = commit_sha[:7]
    return (
        f"Applied the `{stage_name}` curation batch in commit "
        f"[`{short_sha}`]({commit_url}).\n\n"
        "The staged records were applied, committed, and pushed with the "
        "refreshed production-data baseline."
    )


def close_github_issues(
    issue_numbers: list[int], stage_name: str, commit_sha: str
) -> None:
    if shutil.which("gh") is None:
        raise SystemExit(
            "GitHub issue closure requires the GitHub CLI (`gh`) on PATH."
        )

    repo = github_repo_from_remote_url(git_output(["remote", "get-url", "origin"]))
    commit_url = f"https://github.com/{repo}/commit/{commit_sha}"
    comment = issue_close_comment(stage_name, commit_sha, commit_url)
    for issue_number in issue_numbers:
        command = [
            "gh",
            "issue",
            "close",
            str(issue_number),
            "--repo",
            repo,
            "--reason",
            "completed",
            "--comment",
            comment,
        ]
        try:
            subprocess.run(
                command,
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
            )
        except subprocess.CalledProcessError as exc:
            output = (exc.stderr or exc.stdout or "").strip()
            if not output:
                output = (
                    f"`gh issue close` exited with status {exc.returncode} "
                    "without additional output."
                )
            raise SystemExit(
                "Commit/push already completed, but closing "
                f"GitHub issue #{issue_number} failed.\n"
                "Retry manually after resolving the GitHub CLI error:\n"
                f"  {shlex.join(command)}\n\n"
                f"GitHub CLI output:\n{output}"
            ) from exc
        print(f"Closed GitHub issue #{issue_number}.")


def main() -> None:
    args = parse_args()
    manifest_path = manifest_path_from_args(args).resolve()
    manifest = load_manifest(manifest_path)

    if args.commit_message and not manifest.get("applied"):
        raise SystemExit("Cannot commit a staging batch that has not been applied.")
    verify_manifest_hashes(manifest)
    if args.commit_message and not args.skip_verification:
        run_curation_verification()

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
    commit_sha = git_output(["rev-parse", "HEAD"])

    if args.push:
        subprocess.run(
            ["git", "push"],
            cwd=ROOT,
            check=True,
        )
        print("Pushed current branch.")

    if args.close_issue:
        close_github_issues(
            args.close_issue,
            stage_name_from_args(args, manifest_path),
            commit_sha,
        )


if __name__ == "__main__":
    main()
