from abc import ABC, abstractmethod
from neoqm.io import convert

class QMBaseWrapper(ABC):

    def __init__(self):
        self.is_open_shelled = False

    @abstractmethod
    def set_charge(self,charge):
        pass
    
    @abstractmethod            
    def get_energy_and_gradient(self, geometry, point_charges=None):
        # point charge in [n,4]: q,x,y,z
        pass