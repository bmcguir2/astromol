from astromol.database import Database

db = Database()
print(db)
print()

print(f"Loaded {len(db.refs)} references")
print(f"Loaded {len(db.telescopes)} telescopes")
print(f"Loaded {len(db.sources)} sources")
print(f"Loaded {len(db.molecules)} molecules")
print(f"Loaded {len(db.detections)} detections")
print()

for label, mol in sorted(db.molecules.items()):
    print(f"{label}:")
    print(f"  name: {mol.name}")
    print(f"  formula: {mol.formula}")
    print(f"  mass: {mol.mass}")
    print(f"  natoms: {mol.natoms}")
    print(f"  charge: {mol.charge}")
    print(f"  radical: {mol.radical}")
    print(f"  is_linear: {mol.is_linear}")
    if mol.rotcon is not None:
        print(f"  rotcon refs: {[ref.bibcode for ref in mol.rotcon.refs]}")
    if mol.dipole is not None:
        print(f"  dipole refs: {[ref.bibcode for ref in mol.dipole.refs]}")
    print()

print("Detections:")
for detection in db.detections:
    print(f"  {detection}")
