def compute_asymmetry_parameter(A, B, C):
    """
    Returns Ray's asymmetry parameter κ, where:
    κ = (2B - A - C) / (A - C)
    Returns -1 if A and C are both None or degenerate.
    """
    if A is None and C is None:
        return -1.0
    if None in (A, B, C):
        return None
    if A == C:
        return None  # degenerate case
    try:
        return round((2 * B - A - C) / (A - C), 6)
    except ZeroDivisionError:
        return None

def apply_asymmetry(molecule):
    """
    Updates molecule["molecular_properties"]["rays_asymmetry"] using A, B, C.
    """
    props = molecule.get("molecular_properties", {})
    rot = props.get("rotational_constants", {})
    kappa = compute_asymmetry_parameter(rot.get("A"), rot.get("B"), rot.get("C"))
    molecule["molecular_properties"]["rays_asymmetry"] = kappa
    return molecule