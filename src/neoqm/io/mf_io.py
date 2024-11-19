from .mol_io import rdmol2scfmol


def mol2mf_from_chkfile(mol, chkfile, basis="sto-3g"):
    mol_scf = rdmol2scfmol(mol, basis=basis)
    mf = mol_scf.RHF()
    mf.__dict__.update(pyscf.scf.chkfile.load(chkfile, 'scf'))
    return mf
