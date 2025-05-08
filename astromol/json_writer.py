import json
import os

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")

def _get_path(filename):
    return os.path.join(DATA_DIR, filename)

def load_json(filename):
    with open(_get_path(filename), "r", encoding="utf-8") as f:
        return json.load(f)

def save_json(filename, data):
    with open(_get_path(filename), "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)

def append_molecule(entry):
    molecules = load_json("molecules.json")
    existing_ids = {m["id"] for m in molecules}
    if entry["id"] in existing_ids:
        raise ValueError(f"Molecule with ID '{entry['id']}' already exists.")
    molecules.append(entry)
    save_json("molecules.json", molecules)
    print(f"Appended molecule '{entry['id']}' to molecules.json")