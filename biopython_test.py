from Bio.PDB import PPBuilder, calc_angle, calc_dihedral
from Bio.PDB.PDBParser import PDBParser  # https://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ


parser = PDBParser(PERMISSIVE=True)

# Load pdb from benchmark
pdbPath = '../benchmark 5/structures/1A2K_l_b.pdb'
pdbFile = open(pdbPath)
structure = parser.get_structure('1A2K_l_b', pdbFile)

# Iterate over all atoms in a structure
for atom in structure.get_atoms():
    print(atom)

# Iterate over all residues in a model
for residue in structure.get_residues():
    print(residue)

# Calc angle
atoms = structure.get_atoms()
vector1 = next(atoms).get_vector()
vector2 = next(atoms).get_vector()
vector3 = next(atoms).get_vector()
vector4 = next(atoms).get_vector()
angle = calc_angle(vector1, vector2, vector3)
print(angle)

# Calc torsion angles
print(calc_dihedral(vector1, vector2, vector3, vector4))

# Build peptide
ppb = PPBuilder()
for pp in ppb.build_peptides(structure):
    print(pp.get_sequence())
