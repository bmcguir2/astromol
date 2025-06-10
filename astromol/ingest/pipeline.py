from .scaffolder import create_base_molecule
from .composition import apply_elemental_composition
from .asymmetry import apply_asymmetry
from .resolver import apply_name_resolution


def enrich_molecule(molecule):
    """
    Apply all enrichment steps to a molecule JSON dict.
    """
    molecule = apply_elemental_composition(molecule)
    molecule = apply_asymmetry(molecule)
    molecule = apply_name_resolution(molecule)
    return molecule


def create_and_enrich(formula, name=None, state_class=None):
    """
    Creates a new molecule and runs full enrichment.
    """
    mol = create_base_molecule(formula, name=name, state_class=state_class)
    return enrich_molecule(mol)