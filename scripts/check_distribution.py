"""Check the data and curator/test contents of built wheels and source archives."""

from __future__ import annotations

import argparse
from pathlib import Path, PurePosixPath
import tarfile
import zipfile

DATA_FILES = {"molecules.json", "detections.json", "sources.json", "telescopes.json", "references.bib"}
SOURCE_FILES = {
    "scripts/stage_records.py", "scripts/update_data_baseline.py", "scripts/cleanup_stage.py",
    "scripts/check_curation.py", "scripts/generate_outputs.py", "scripts/check_distribution.py",
    "tests/baselines/production_data.json", "tests/baseline.py", "tests/rendering_assertions.py",
    "curation/README.md", "curation/templates/molecule.yaml", "curation/templates/detection.yaml",
    "curation/templates/source.yaml", "curation/templates/telescope.yaml", "docs/conf.py", "SPEC.md",
}


def check_distribution(path: Path) -> None:
    if path.suffix == ".whl":
        with zipfile.ZipFile(path) as archive:
            names = {name for name in archive.namelist() if not name.endswith("/")}
    else:
        with tarfile.open(path) as archive:
            names = {PurePosixPath(member.name).relative_to(PurePosixPath(member.name).parts[0]).as_posix()
                     for member in archive.getmembers() if member.isfile()}
        missing = SOURCE_FILES - names
        if missing:
            raise ValueError(f"{path}: missing source workflow files: {sorted(missing)}")
    packaged_data = {PurePosixPath(name).name for name in names if name.startswith("astromol/data/")}
    if packaged_data != DATA_FILES:
        raise ValueError(f"{path}: expected only canonical production data, found {sorted(packaged_data)}")
    forbidden = {name for name in names if name.startswith("curation/staging/") or ".preview." in name
                 or "_stage_report." in name or "_stage_manifest." in name or ".apply-backup/" in name}
    if forbidden:
        raise ValueError(f"{path}: temporary curation artifacts included: {sorted(forbidden)}")
    print(f"Checked {path.name}: canonical data and complete source workflow.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archives", type=Path, nargs="+")
    for path in parser.parse_args().archives:
        check_distribution(path)


if __name__ == "__main__":
    main()
