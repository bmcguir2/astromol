import requests
from rdkit import Chem
import selfies


def fetch_smiles_from_pubchem(name):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/CanonicalSMILES/JSON"
    try:
        r = requests.get(url, timeout=10)
        r.raise_for_status()
        return r.json()["PropertyTable"]["Properties"][0]["CanonicalSMILES"]
    except Exception as e:
        raise ValueError(f"PubChem lookup failed for '{name}': {e}")

def apply_name_resolution(molecule):
    ids = molecule.get("identifiers", {})
    name = ids.get("name")

    if not ids.get("canonical_smiles") and name:
        try:
            ids["canonical_smiles"] = fetch_smiles_from_pubchem(name)
        except Exception as e:
            molecule.setdefault("notes", "")
            molecule["notes"] += f" [resolver warning: {e}]"
            return molecule

    smiles = ids.get("canonical_smiles")
    if not smiles:
        molecule.setdefault("notes", "")
        molecule["notes"] += " [resolver warning: no SMILES available]"
        return molecule

    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            raise ValueError("Invalid SMILES string")

        inchi = Chem.MolToInchi(mol)
        inchikey = Chem.InchiToInchiKey(inchi)
        selfies_str = selfies.encoder(smiles)

        ids["selfies"] = selfies_str
        ids["inchi"] = inchi
        ids["inchikey"] = inchikey

    except Exception as e:
        molecule.setdefault("notes", "")
        molecule["notes"] += f" [resolver error: {e}]"

    return molecule