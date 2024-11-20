from operator import ge
import os
import time
import copy
import shutil
import numpy as np
from neoqm import unit
from neoqm.utils import run_command, RuntimeTempDir
from neoqm.base import QMBaseWrapper

class OrcaWrapper(QMBaseWrapper):

    def __init__(self,
                 method='RKS',
                 xc='b3lyp',
                 basis='6-31g*',
                 directory=None,
                 **kwargs):

        self.orca_exe = kwargs.get("orca", shutil.which("orca"))
        self.orca_job_name =kwargs.get( "orca_job_name","run")
        #  self.file_dic = {}
        self.method = method
        self.xc = xc
        self.basis = basis
        self.directory = directory
        self.threads = kwargs.get("threads", 10)
        self.charge = kwargs.get("charge", 0)
        self.multiplicity = kwargs.get("multiplicity", 1)
        
    def set_charge(self,charge):
        self.charge = charge
    @property
    def engine_name(self):
        return "orca"
    def get_command(self,directory):
        return [self.orca_exe, os.path.join(directory,"run.in")]

    def edit_template(self, geom, first_line_add):
        line1_ls = [self.method, self.xc, self.basis]
        line1_ls.remove(None)
        if first_line_add is not None:
            if isinstance(first_line_add, list):
                for addition in first_line_add:
                    line1_ls.append(addition)
            else:
                line1_ls.append(addition)
        _line1 = ' '.join(line1_ls)
        _line1 = '! ' + _line1
        template = [
            _line1, "%PAL NPROCS {} END".format(self.threads),
            "* XYZ {} {}".format(self.charge, self.multiplicity)
        ]
        template = "\n".join(template)
        template += "\n"
        template += geom
        template += "*"
        return template

    def execute(self, directory):
        start_time = time.time()

        stdout_file = os.path.join(directory,
                                   "{}.out".format(self.orca_job_name))
        stderr_file = os.path.join(directory,
                                   "{}.err".format(self.orca_job_name))

        runtime = run_command(self.command,
                              stdout_file=stdout_file,
                              stderr_file=stderr_file)
        time_used = time.time() - start_time

        if runtime:
            raise RuntimeError("Orca Runtime Error")
        with open(stdout_file, 'w') as f:
            f.write(runtime.stdout.decode("utf-8"))
        with open(stderr_file, 'w') as f:
            f.write(runtime.stderr.decode("utf-8"))
        return True

    def run(self, geom, directory, first_line_add=None):
        template = self.edit_template(geom, first_line_add)
        print('working dir for orca job:{}\n'.format(directory))
        with open(os.path.join(directory, "{}.in".format(self.orca_job_name)),
                  "w") as f:
            f.write(template)

        self.execute(directory)

    @staticmethod
    def parse_energy(filename):
        out = []
        with open(filename, 'r') as f:
            lines = f.readlines()
        for line in lines:
            if 'FINAL SINGLE POINT ENERGY' in line:
                out.append(line.split()[-1])
        return out

    @staticmethod
    def parse_gradient(filename, atom_nums=None):
        with open(filename, 'r') as f:
            lines = f.readlines()
        for i in range(len(lines)):
            if '# The current gradient' in lines[i]:
                break
        line0 = i + 2
        line1 = i + 2 + atom_nums * 3

        grad = np.array([float(x.strip(' \n')) for x in lines[line0:line1]])
        grad.reshape(atom_nums, 3)
        return grad
    

    def get_energy_and_gradient(self, geometry, point_charges=None,unit_in="atomic"):
        # point charge in [n,4]: q,x,y,z
        # template with {geometry} \n Q q x y z
        # geom should be xyz content string
        result = {}
        geometry = copy.deepcopy(geometry)
        # ensure the last line is a newline
        geometry = geometry.strip('\n') + '\n'
        with RuntimeTempDir(self.directory) as tmpdir:
            if point_charges is not None:
                for point_charge in point_charges:
                    geometry += "Q {} {} {} {}\n".format(*point_charge)
            self.run(geometry, tmpdir.name, first_line_add=['ENGRAD'])
            energy_path = os.path.join(tmpdir.name,
                                       '{}.out'.format(self.orca_job_name))
            grad_path = os.path.join(tmpdir.name,
                                     '{}.engrad'.format(self.orca_job_name))
            energy = self.parse_energy(energy_path)
            gradient = self.parse_gradient(
                grad_path, atom_nums=len(geometry.strip('\n').split('\n')))
        energy = energy * unit.hartree * unit.avogadro_constant
        gradient = gradient * unit.hartree * unit.avogadro_constant / unit.bohr
        if unit_in in {"openmm","ttk"}:
            energy = energy.to(unit.kilojoule / unit.mole)
            gradient = gradient.to(unit.kilojoule / unit.mole / unit.nanometer)
        return energy ,gradient
