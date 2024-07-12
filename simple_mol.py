
from rdkit import Chem
import numpy as np


class SimpleMol:
    """
        A class for storing molecules in basic python datatypes
    
    """

    def __init__(self, smiles: str = None, rdmol: Chem.Mol = None) -> None:
        """
            Initialise a SimpleMol using a smiles string, or rdkit Chem.Mol object as input
        """

        self.elements = []
        self.bonds = []

        if smiles is not None:
            rdmol = Chem.MolFromSmiles(smiles)

        elements, bonds = rdmol_to_simple(rdmol)
        self.elements = elements
        self.bonds = bonds

    def to_rdmol(self):

        # create empty editable mol object
        mol = Chem.RWMol()

        # add atoms to mol and keep track of index
        node_to_idx = {}
        for i in range(len(self.elements)):
            a = Chem.Atom(self.elements[i])
            molIdx = mol.AddAtom(a)
            node_to_idx[i] = molIdx

        # add bonds between adjacent atoms
        for ix, row in enumerate(self.bonds):
            for iy, bond in enumerate(row):

                # only traverse half the matrix
                if iy <= ix:
                    continue

                # add relevant bond type (there are many more of these)
                if bond == 0:
                    continue
                elif bond == 1:
                    bond_type = Chem.rdchem.BondType.SINGLE
                    mol.AddBond(node_to_idx[ix], node_to_idx[iy], bond_type)
                elif bond == 2:
                    bond_type = Chem.rdchem.BondType.DOUBLE
                    mol.AddBond(node_to_idx[ix], node_to_idx[iy], bond_type)

        # Convert RWMol to Mol object
        mol = mol.GetMol()     
        return mol

    def to_smiles(self) -> str:
        """ Output simple mol as smiles string """

        mol = self.to_rdmol()
        return Chem.MolToSmiles(mol)
        

def read_sdf_file(file: str):
    """
        Reads in a series of smiles strings from a file, converting each smiles
            string to a SimpleMol object
    """

    mols = []
    with Chem.SDMolSupplier(file) as suppl:
        for mol in suppl:
            mols.append(mol)

    simplemols = []
    for mol in mols:
        simplemol = SimpleMol(rdmol=mol)
        simplemols.append(simplemol)

    return simplemols

def read_smiles_file(file: str):
    """
        Reads in a series of smiles strings from a file, converting each smiles
            string to a SimpleMol object
    """

    smiles = []
    with open(file, 'r') as f:
        for line in f:
            smiles.append(line.strip())

    simplemols = []
    for smile in smiles:
        simplemol = SimpleMol(smiles=smile)
        simplemols.append(simplemol)

    return simplemols




def rdmol_to_simple(rdmol):
    """
        take basic molecular information out of an rdkit molecule
    """

    num_atoms = rdmol.GetNumAtoms()
    bonds = np.zeros((num_atoms, num_atoms))
    elements = []

    for atom in rdmol.GetAtoms():
        element = atom.GetSymbol()
        elements.append(element)

    for bond in rdmol.GetBonds():
        index1 = bond.GetBeginAtomIdx()
        index2 = bond.GetEndAtomIdx()
        bonds[index1][index2] = bond.GetBondTypeAsDouble()

    return elements, bonds