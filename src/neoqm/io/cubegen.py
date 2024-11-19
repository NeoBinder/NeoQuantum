import numpy as np
from pyscf import lib, gto, df
from pyscf.dft import numint
from pyscf.pbc.gto import Cell
from pyscf.tools.cubegen import Cube


def generate_cube_specifiedbox(mol, resolution, origin, margin, extent):
    return Cube(mol,
                resolution=resolution,
                origin=origin,
                margin=margin,
                extent=extent)


def generate_cube(mol):
    return Cube(mol)


def density(cube, dm):
    """Calculates electron density and write out in cube format.
    Args:
        cube : Cube
            Cube to calculate the eletron density for
        dm : ndarray
            Density matrix of molecule.
    """
    mol = cube.mol
    GTOval = 'GTOval'
    if isinstance(mol, Cell):
        GTOval = 'PBC' + GTOval

    # Compute density on the .cube grid
    coords = cube.get_coords()
    ngrids = cube.get_ngrids()
    blksize = min(8000, ngrids)
    rho = np.empty(ngrids)
    for ip0, ip1 in lib.prange(0, ngrids, blksize):
        ao = mol.eval_gto(GTOval, coords[ip0:ip1])
        rho[ip0:ip1] = numint.eval_rho(mol, ao, dm)
    rho = rho.reshape(cube.nx, cube.ny, cube.nz)
    #  # Write out density to the .cube file
    #  cc.write(rho, outfile, comment='Electron density in real space (e/Bohr^3)')
    return rho


def mep(cube, dm):
    """Calculates the molecular electrostatic potential (MEP) and write out in
    cube format.
    Args:
        cube : Cube
            cube to calculate the electron density for.
        dm : ndarray
            Density matrix of molecule.
    """

    mol = cube.mol
    coords = cube.get_coords()

    # Nuclear potential at given points
    Vnuc = 0
    for i in range(mol.natm):
        r = mol.atom_coord(i)
        Z = mol.atom_charge(i)
        rp = r - coords
        Vnuc += Z / np.einsum('xi,xi->x', rp, rp)**.5

    # Potential of electron density
    Vele = np.empty_like(Vnuc)
    for p0, p1 in lib.prange(0, Vele.size, 600):
        fakemol = gto.fakemol_for_charges(coords[p0:p1])
        ints = df.incore.aux_e2(mol, fakemol)
        Vele[p0:p1] = np.einsum('ijp,ij->p', ints, dm)

    MEP = Vnuc - Vele  # MEP at each point
    MEP = MEP.reshape(cube.nx, cube.ny, cube.nz)

    # Write the potential
    #  cc.write(MEP, outfile, 'Molecular electrostatic potential in real space')
    return MEP
