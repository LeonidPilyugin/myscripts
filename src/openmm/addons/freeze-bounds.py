from addon import OpenmmAddOn, OpenmmSimulation
from gi.repository import AmlCore, AmlBasicTypes, AmlParticles
from pathlib import Path
import os
import openmm
import numpy as np

class FreezeBoundsActionParams(AmlCore.ActionParams):
    def __init__(self):
        super().__init__()
        self.top = 0
        self.bottom = 0

    def do_copy(self):
        res = FreezeBoundsActionParams()
        res.top = self.top
        res.bottom = self.bottom
        return res

class FreezeBoundsAction(AmlCore.Action):
    def __init__(self):
        super().__init__()

    def do_get_params_error_message(self, params):
        return ""

    def do_perform(self, data : AmlCore.DataCollection):
        logger = data.get_element("logger").logger
        params = self.get_params()

        particles = data.get_element("repr.particles")
        openmm_obj = data.get_element("openmm")

        system = openmm_obj.context.getSystem()
        state = openmm_obj.context.getState(getPositions=True, getVelocities=True, enforcePeriodicBox=False)
        n = system.getNumParticles()
        positions = state.getPositions(asNumpy=True)
        velocities = state.getVelocities(asNumpy=True)
        masses = np.asarray(particles.get_prop("mass").get_arr())

        logger.info("Getting new masses")

        for i in range(n):
            if positions[i][2] > params.top * openmm_obj.unit.length or positions[i][2] < params.bottom * openmm_obj.unit.length:
                system.setParticleMass(i, 0.0)
                masses[i] = 0.0

        logger.info("Reinitializing OpenMM")

        openmm_obj.context.reinitialize()
        openmm_obj.context.setPositions(positions)
        openmm_obj.context.setVelocities(velocities)

        logger.info("Updating repr DataObject")

        particles.set_prop("mass", AmlParticles.Float64PerParticleProperty.from_array(masses.tolist()))


class FreezeBounds(OpenmmAddOn):
    def __init__(self):
        self.action = FreezeBoundsAction()
        self.params = FreezeBoundsActionParams()

    def init_dict(self, dictionary):
        self.params.top = dictionary["top"]
        self.params.bottom = dictionary["bottom"]
        self.action.set_params(self.params)

    def get_action(self):
        return self.action
