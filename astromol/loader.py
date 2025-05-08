import json
import os

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE_DIR, "..", "data")


def _load_json(filename):
    path = os.path.join(DATA_DIR, filename)
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def load_molecules():
    return _load_json("molecules.json")


def load_sources():
    return _load_json("sources.json")


def load_telescopes():
    return _load_json("telescopes.json")


def load_references():
    return _load_json("references.json")