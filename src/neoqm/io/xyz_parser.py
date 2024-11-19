import numpy as np
from mendeleev import element as mendeleev_element
from mendeleev import Element

symbol_element_dict = {}

__all__ = ["XYZFile"]


def element_from_symbol(symbol):
    if symbol not in symbol_element_dict:
        element = mendeleev_element(symbol)
        symbol_element_dict[symbol] = element
    return symbol_element_dict[symbol]


class XYZFile():

    def __init__(self, atoms, coords):
        # in nm
        self.atoms = atoms
        self.coords = coords

    @classmethod
    def from_file(cls, filepath):
        with open(filepath, "r") as f:
            count = int(f.readline().strip())
            _ = f.readline()
            atoms = []
            coords = []
            data = f.readlines()
            data = map(lambda x: x.strip().split(), data)
            for line in data:
                atoms.append(element_from_symbol(line[0]))
                coords.append( line[1:])
            return cls(atoms, np.vstack(coords).astype(np.float64))
