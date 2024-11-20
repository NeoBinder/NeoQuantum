import numpy as np
from pyscf import gto, scf, grad
from neoqm import unit
from neoqm.base import QMBaseWrapper



class PyscfWrapper(QMBaseWrapper):
    def __init__(
        self,
        method="scf",
        charge=0.0,
        multiplicity=1,
        reference="rhf",
        basis="STO-3G",
        **kwargs,
    ):

        super().__init__()
        self.basis = basis.lower()
        self.set_charge(charge)

        self.method = method
        self.multiplicity = multiplicity
        #  self.external_charges = None

        self.qm_param = kwargs
        self.qm_param["reference"] = reference
        self.qm_param["FAIL_ON_MAXITER"] = False

    @property
    def engine_name(self):
        return "pyscf"

    def set_charge(self, charge, quiet=True):
        if hasattr(self, "mol"):
            self.mol.charge = charge
        else:
            self.mol = gto.Mole()
            self.mol.basis = self.basis
            self.mol.unit = "Angstrom"
            self.mol.max_memory = 1e4  # 10G
            self.mol.charge = charge
        if quiet:
            self.mol.verbose = 0
        self.mol.build()
        self.scanner = scf.RHF(self.mol).apply(grad.RHF).as_scanner()

    def get_energy_and_gradient(self, geometry,unit_in="atomic"):
        """
        Calculate the energy and gradient of the molecule at the given geometry.

        Parameters:
        geometry (numpy.ndarray): The molecular geometry in Angstrom.

        Returns:
        dict: A dictionary containing the energy and gradients.
        """
        energy, gradient = self.scanner(geometry)        
        energy = energy * unit.hartree * unit.avogadro_constant
        gradient = gradient * unit.hartree * unit.avogadro_constant / unit.bohr
        if unit_in in {"openmm","ttk"}:
            energy = energy.to(unit.kilojoule / unit.mole)
            gradient = gradient.to(unit.kilojoule / unit.mole / unit.nanometer)
        return energy , gradient
