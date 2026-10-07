"""Small, explicit provenance records for generated scientific products."""

from __future__ import annotations

import hashlib
from importlib import metadata
from pathlib import Path
import platform
import subprocess

DATA_FILES = ("molecules.json", "detections.json", "sources.json", "telescopes.json", "references.bib")
DEPENDENCIES = ("astromol", "bibtexparser", "matplotlib", "molmass", "numpy", "python-pptx", "scipy")


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_revision() -> dict:
    root = Path(__file__).resolve().parents[1]
    if not (root / ".git").exists():
        return {"commit": None, "dirty": None}
    try:
        def git(*args):
            return subprocess.run(["git", "-C", str(root), *args], check=True,
                                  capture_output=True, text=True).stdout.strip()
        return {"commit": git("rev-parse", "HEAD"), "dirty": bool(git("status", "--porcelain"))}
    except (OSError, subprocess.SubprocessError):
        return {"commit": None, "dirty": None}


def environment_versions() -> dict:
    versions = {"python": platform.python_version()}
    for name in DEPENDENCIES:
        try:
            versions[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            versions[name] = None
    return versions


def generation_provenance(data_dir: Path) -> dict:
    package = Path(__file__).resolve().parent
    return {
        "source": source_revision(),
        "data_sha256": {name: digest(data_dir / name) for name in DATA_FILES},
        "code_sha256": {p.relative_to(package).as_posix(): digest(p) for p in sorted(package.rglob("*.py"))},
        "versions": environment_versions(),
        "platform": platform.platform(),
    }
