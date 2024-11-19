import os
from qmkit.qmcalc import mol2mf
from qmkit.io import file2mol, cubegen
from qmkit import utils


def mol2grid(fpath, density=True, mep=True):
    mol = file2mol(fpath)
    mol_scf, mf = mol2mf(mol)

    cube = cubegen.generate_cube(mol_scf)
    output = []
    if density:
        output.append(cubegen.density(cube, mf.make_rdm1()))
    if mep:
        output.append(cubegen.mep(cube, mf.make_rdm1()))
    return output


def test_mol2grid():
    ligand_fpath = os.path.join(utils.DATADIR, "1a30", 'ligand.mol2')
    print(mol2grid(ligand_fpath))
