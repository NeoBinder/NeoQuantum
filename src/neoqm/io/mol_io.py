import os
from rdkit import Chem
import pyscf

__all__ = ["file2mol", "rdmol2xyz", "rdmol2scfmol"]


def file2mol(filepath):
    filename, file_extension = os.path.splitext(filepath)
    if file_extension not in (".pdb", ".mol2", ".mol", ".sdf"):
        raise NotImplementedError('Unrecognized file format with \"' +
                                  file_extension + '\"!!!\n')

    if file_extension == '.pdb':
        mol = Chem.MolFromPDBFile(filepath, removeHs=False)
    elif file_extension == '.mol2':
        mol = Chem.MolFromMol2File(filepath, removeHs=False)
    elif file_extension == '.mol':
        mol = Chem.MolFromMolFile(filepath, removeHs=False)
    elif file_extension == '.sdf':
        supp = Chem.ForwardSDMolSupplier(filepath, removeHs=False)
        mol = next(supp)
    return mol


def rdmol2xyz(mol_rdkit):
    charge = Chem.GetFormalCharge(mol_rdkit)
    pos_list = mol_rdkit.GetConformer().GetPositions()
    atom_list = [a.GetSymbol() for a in mol_rdkit.GetAtoms()]
    mol_str = ''
    for i in range(len(atom_list)):
        atom = atom_list[i]
        x = pos_list[i][0]
        y = pos_list[i][1]
        z = pos_list[i][2]
        mol_str += atom + ' ' + str(x) + ' ' + str(y) + ' ' + str(z) + "\n"
    return mol_str


def rdmol2scfmol(mol_rdkit, basis="sto-3g"):
    mol_xyz = rdmol2xyz(mol_rdkit)
    charge = Chem.GetFormalCharge(mol_rdkit)

    mol_scf = pyscf.M(atom=mol_xyz,
                      unit='Angstrom',
                      basis=basis,
                      charge=charge)
    return mol_scf
