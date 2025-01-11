from addon import OpenmmAddOn, OpenmmSimulation
from gi.repository import AmlCore, AmlBasicTypes, AmlParticles
from pathlib import Path
import os
import openmm
import numpy as np

class DisplaceBoundsActionParams(AmlCore.ActionParams):
    def __init__(self):
        super().__init__()
        self.top = 0
        self.bottom = 0
        self.magnitude = 0

    def do_copy(self):
        res = DisplaceBoundsActionParams()
        res.top = self.top
        res.bottom = self.bottom
        res.magnitude = self.magnitude
        return res

class DisplaceBoundsAction(AmlCore.Action):
    def __init__(self):
        super().__init__()

    def do_get_params_error_message(self, params):
        return ""

    def do_perform(self, data : AmlCore.DataCollection):
        params = self.get_params()

        particles = data.get_element("repr.particles")
        openmm_obj = data.get_element("openmm")

        system = openmm_obj.context.getSystem()
        state = openmm_obj.context.getState(getPositions=True, getVelocities=True, enforcePeriodicBox=False)
        n = system.getNumParticles()
        positions = state.getPositions(asNumpy=True)
        velocities = state.getVelocities(asNumpy=True)

        if params.top is None:
            assert params.bottom is not None
            for i in range(n):
                if positions[i][2] < params.bottom * openmm_obj.unit.length:
                    positions[i][0] += params.magnitude * openmm_obj.unit.length
        else:
            assert params.bottom is None
            for i in range(n):
                if positions[i][2] > params.top * openmm_obj.unit.length:
                    positions[i][0] += params.magnitude * openmm_obj.unit.length


        openmm_obj.context.setPositions(positions)
        openmm_obj.context.setVelocities(velocities)

        particles.set_prop("x", AmlParticles.Float64PerParticleProperty.from_array(positions[:,0].tolist()))


class DisplaceBounds(OpenmmAddOn):
    def __init__(self):
        self.action = DisplaceBoundsAction()
        self.params = DisplaceBoundsActionParams()

    def init_dict(self, dictionary):
        self.params.top = dictionary.get("top", None)
        self.params.bottom = dictionary.get("bottom", None)
        self.params.magnitude = dictionary["magnitude"]
        self.action.set_params(self.params)

    def get_action(self):
        return self.action
