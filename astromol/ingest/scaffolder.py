from datetime import datetime
import re

def create_base_molecule(formula, name=None, state_class=None):
    """
    Create a minimally valid molecule JSON structure with properly typed defaults.
    Derives ID from formula, but prepends positional isomer prefix if present in name.
    """
    base_id = formula
    if name:
        match = re.match(r"^(\d+|cis|trans|iso|sec|tert|n|t|l|d)-", name, re.IGNORECASE)
        if match:
            base_id = f"{match.group(0)}{formula}"

    return {
        "id": base_id,
        "state_class": state_class or ["canonical_form"],
        "identifiers": {
            "name": name or formula,
            "formula": formula,
            "canonical_smiles": "",
            "selfies": "",
            "inchi": "",
            "inchikey": "",
            "other_names": []
        },
        "tags": [],
        "latex": {
            "table_formula": base_id,
            "label": f"mol:{base_id}"
        },
        "detections": [],
        "manuscript_text": {"latex": ""},
        "molecular_properties": {
            "natoms": None,
            "mass": None,
            "rotational_constants": {
                "A": None,
                "B": None,
                "C": None
            },
            "dipole_moment": {
                "a": None,
                "b": None,
                "c": None
            },
            "rays_asymmetry": None,
            "element_counts": {}
        },
        "notes": "",
        "history": [
            {
                "date": datetime.utcnow().isoformat(timespec="seconds") + "Z",
                "version": "0.1.0",
                "user": "bmcguire",
                "note": f"Initial creation of {base_id}"
            }
        ]
    }