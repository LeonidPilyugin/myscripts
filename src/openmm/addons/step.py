from addon import OpenmmAddOn, OpenmmSimulation
from gi.repository import AmlCore, AmlParticles, AmlBasicTypes
import numpy as np
import openmm
import openmm.unit
import scipy
import scipy.constants

class PerformStepActionParams(AmlCore.ActionParams):
    def __init__(self):
        super().__init__()
        self.number = 0
        self.mean = False

    def do_copy(self):
        res = PerformStepActionParams()
        res.number = self.number
        res.mean = self.mean
        return res

class PerformStepAction(AmlCore.Action):
    def __init__(self):
        super().__init__()

    def do_get_params_error_message(self, params):
        if not isinstance(params, PerformStepActionParams):
            return "Invalid params type"
        if not isinstance(params.number, int) or params.number < 0:
            return "Invalid number"
        return ""

    def do_perform(self, data : AmlCore.DataCollection):
        params = self.get_params()
        openmm_object = data.get_element("openmm")

        state = openmm_object.context.getState(
            {
                "getPositions": True,
            }
        )

        u = t = T = 0
        positions = np.zeros_like(state.getPositions(asNumpy=True))
        velocities = np.zeros_like(positions)
        masses = np.asarray(data.get_element("repr.particles").get_prop("mass").get_arr())
        masses /= scipy.constants.N_A * 1000

        if (params.mean):
            for _ in range(params.number):
                openmm_object.integrator.step(1)
                state = openmm_object.context.getState(
                    getPositions=True,
                    getVelocities=True,
                    enforcePeriodicBox=False,
                    getForces=True,
                    getEnergy=True,
                )
                p = state.getPositions(asNumpy=True)
                v = state.getVelocities(asNumpy=True)
                positions = np.add(positions, p.value_in_unit(openmm_object.unit.length))
                velocities = np.add(velocities, v.value_in_unit(openmm_object.unit.velocity))
                u += state.getPotentialEnergy().value_in_unit(openmm_object.unit.energy)
                t += state.getKineticEnergy().value_in_unit(openmm_object.unit.energy)
                T += (masses * np.sum(v.value_in_unit(openmm.unit.meter / openmm.unit.second) ** 2, axis=1)).sum() / scipy.constants.k / np.count_nonzero(masses) / 3

            positions /= params.number
            velocities /= params.number
            u /= params.number
            t /= params.number
            T /= params.number
        else:
            openmm_object.integrator.step(params.number)
            state = openmm_object.context.getState(
            getPositions=True,
                getVelocities=True,
                enforcePeriodicBox=False,
                getForces=True,
                getEnergy=True,
            )

            positions = state.getPositions(asNumpy=True)
            velocities = state.getVelocities(asNumpy=True)
            u = state.getPotentialEnergy().value_in_unit(openmm_object.unit.energy)
            t = state.getKineticEnergy().value_in_unit(openmm_object.unit.energy)
            T = (masses * np.sum(velocities.value_in_unit(openmm.unit.meter / openmm.unit.second) ** 2, axis=1)).sum() / scipy.constants.k / np.count_nonzero(masses) / 3
            positions = positions.value_in_unit(openmm_object.unit.length)
            velocities = velocities.value_in_unit(openmm_object.unit.velocity)

        data.get_element("repr.thermo.temperature").set_val(T)
        data.get_element("repr.thermo.potential").set_val(u)
        data.get_element("repr.thermo.kinetic").set_val(t)
        data.get_element("repr.thermo.total").set_val(u + t)

        step = data.get_element("repr.step").get_val()
        data.get_element("repr.step").set_val(step + params.number)

        particles = data.get_element("repr.particles")
        particles.set_prop("x", AmlParticles.Float64PerParticleProperty.from_array(positions[:,0].tolist()))
        particles.set_prop("y", AmlParticles.Float64PerParticleProperty.from_array(positions[:,1].tolist()))
        particles.set_prop("z", AmlParticles.Float64PerParticleProperty.from_array(positions[:,2].tolist()))
        particles.set_prop("vx", AmlParticles.Float64PerParticleProperty.from_array(velocities[:,0].tolist()))
        particles.set_prop("vy", AmlParticles.Float64PerParticleProperty.from_array(velocities[:,1].tolist()))
        particles.set_prop("vz", AmlParticles.Float64PerParticleProperty.from_array(velocities[:,2].tolist()))


class PerformStep(OpenmmAddOn):
    def __init__(self):
        self.action = PerformStepAction()
        self.params = PerformStepActionParams()

    def init_dict(self, dictionary):
        self.params.number = dictionary["number"]
        if "mean" in dictionary: self.params.mean = dictionary["mean"]
        
        self.action.set_params(self.params)

    def get_action(self):
        return self.action


