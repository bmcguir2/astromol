import json

import pytest

import astromol.database as database_module


def write_minimal_data_dir(path, *, telescopes=None, sources=None):
    """Write the smallest data directory needed to instantiate Database."""
    (path / "references.bib").write_text("")
    (path / "telescopes.json").write_text(json.dumps(telescopes or []))
    (path / "sources.json").write_text(json.dumps(sources or []))
    (path / "molecules.json").write_text("[]")
    (path / "detections.json").write_text("[]")


def test_database_rejects_duplicate_telescope_nicks(tmp_path, monkeypatch):
    write_minimal_data_dir(
        tmp_path,
        telescopes=[
            {
                "name": "Example Telescope One",
                "nick": "Example",
                "shortname": "Example One",
                "type": "Single Dish",
                "wavelength": ["cm"],
            },
            {
                "name": "Example Telescope Two",
                "nick": "Example",
                "shortname": "Example Two",
                "type": "Single Dish",
                "wavelength": ["mm"],
            },
        ],
    )
    monkeypatch.setattr(database_module, "DATA_DIR", tmp_path)

    with pytest.raises(ValueError, match="Duplicate telescope nick: Example"):
        database_module.Database()


def test_database_rejects_duplicate_source_nicks(tmp_path, monkeypatch):
    write_minimal_data_dir(
        tmp_path,
        sources=[
            {
                "name": "Example Source One",
                "nick": "ExampleSource",
                "type": "Other",
            },
            {
                "name": "Example Source Two",
                "nick": "ExampleSource",
                "type": "Other",
            },
        ],
    )
    monkeypatch.setattr(database_module, "DATA_DIR", tmp_path)

    with pytest.raises(ValueError, match="Duplicate source nick: ExampleSource"):
        database_module.Database()
