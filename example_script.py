
from simple_mol import SimpleMol, read_smiles_file, read_sdf_file


## Example script to read in molecules from file

file = "molecules/basic_mols.smi"

simple_mols = read_smiles_file(file)

for mol in simple_mols:
    print(mol.elements)
    print(mol.bonds)
    print(mol.to_smiles())
