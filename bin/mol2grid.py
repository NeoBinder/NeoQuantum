from qmkit.qmcalc import mol2mf
from qmkit.io import file2mol, cubegen



def mol2grid(fpath, density=True, mep=True):
    if density_outpath is None and mep_outpath is None:
        raise RuntimeError("please specify output")
    mol = file2mol(fpath)
    mol_scf, mf = mol2mf(mol)

    cube = cubegen.generate_cube(mol_scf)

    output = []
    if density:
        output.append(cubegen.density(cube, mf.make_rdm1()))
    if mep:
        output.append(cubegen.mep(cube, mf.make_rdm1()))
    return output

