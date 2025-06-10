from periodictable import formula as pt_formula


def apply_elemental_composition(molecule):
    """
    Uses periodictable.formula() to compute element_counts, natoms, and mass.
    """
    formula_str = molecule.get("identifiers", {}).get("formula")
    if not formula_str:
        return molecule

    try:
        f = pt_formula(formula_str)
    except Exception as e:
        molecule["molecular_properties"]["element_counts"] = {}
        molecule["molecular_properties"]["natoms"] = None
        molecule["molecular_properties"]["mass"] = None
        molecule.setdefault("notes", "")
        molecule["notes"] += f" [periodictable error: {e}]"
        return molecule

    element_counts = {str(e): int(n) for e, n in f.atoms.items()}
    molecule["molecular_properties"]["element_counts"] = element_counts
    molecule["molecular_properties"]["natoms"] = sum(element_counts.values())
    molecule["molecular_properties"]["mass"] = f.mass
    return molecule
