import pyscf

from neoqm.io import rdmol2scfmol
from pyscf.geomopt.geometric_solver import optimize


def mol2mf(mol,
           solv=None,
           solv_eps=78.3553,
           basis="sto-3g",
           opt_cycles=0,
           scf_cycles=50):
    mol_scf = rdmol2scfmol(mol, basis=basis)

    if solv is None:
        mf = mol_scf.RHF(max_cycle=scf_cycles)
    else:
        mf = mol_scf.RHF(max_cycle=scf_cycles).ddCOSMO()
        mf.with_solvent.eps = solv_eps

    opt = (opt_cycles != 0)

    if opt:
        mol_eq = optimize(mf, maxsteps=opt_cycles)
        result = mol_eq.RHF().run()
    else:
        result = mf.run()
    return mol_scf, result
