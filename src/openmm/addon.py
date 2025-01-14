from gi.repository import AmlCore
from abc import abstractmethod, ABCMeta
import logging
import openmm
import openmm.unit

class Unit:
    def __init__(self):
        self.length = openmm.unit.angstrom
        self.time = openmm.unit.picosecond
        self.velocity = self.length / self.time
        self.energy = openmm.unit.kilojoule_per_mole
        self.mass = openmm.unit.atom_mass_units


class Logger(AmlCore.DataObject):
    def __init__(self):
        super().__init__()
        self.logger = None

    def do_copy(self):
        res = Logger()
        res.logger = self.logger
        return res


class OpenmmSimulation(AmlCore.DataObject):
    def __init__(self):
        super().__init__()
        self.context = None
        self.integrator = None
        self.unit = Unit()

    def do_copy(self):
        res = OpenmmSimulation()
        res.context = self.context
        res.integrator = self.integrator
        return res

class OpenmmAddOn(metaclass=ABCMeta):
    @abstractmethod
    def init_dict(self, dictionary: dict):
        pass

    @abstractmethod
    def get_action(self) -> AmlCore.Action:
        pass
