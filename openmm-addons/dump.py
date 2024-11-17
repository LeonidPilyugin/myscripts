from gi.repository import Aml
from openmm import unit
from pathlib import Path
import os

IO = Aml.Io.create(Aml.LammpsTextDumpParser())

def main(step, simulation, data, *args, **kwargs):
    simulation.root.joinpath("nonmean").mkdir(parents=True, exist_ok=True)

    frame = simulation.frame.copy()

    system = simulation.context.getSystem()
    state = simulation.get_state()
    n = system.getNumParticles()
    positions = state.getPositions(asNumpy=True)
    velocities = state.getVelocities(asNumpy=True)

    x = positions[:,0].tolist()
    y = positions[:,1].tolist()
    z = positions[:,2].tolist()
    vx = velocities[:,0].tolist()
    vy = velocities[:,1].tolist()
    vz = velocities[:,2].tolist()
    types = frame.atoms.get_prop("type").get_arr()
    masses = frame.atoms.get_prop("mass").get_arr()
    ids = list(range(n))

    self.frame.atoms.set_prop("x", Aml.DoublePerAtomProperty.from_array(x))
    self.frame.atoms.set_prop("y", Aml.DoublePerAtomProperty.from_array(y))
    self.frame.atoms.set_prop("z", Aml.DoublePerAtomProperty.from_array(z))
    self.frame.atoms.set_prop("vx", Aml.DoublePerAtomProperty.from_array(vx))
    self.frame.atoms.set_prop("vy", Aml.DoublePerAtomProperty.from_array(vy))
    self.frame.atoms.set_prop("vz", Aml.DoublePerAtomProperty.from_array(vz))
    self.frame.atoms.set_prop("mass", Aml.DoublePerAtomProperty.from_array(masses))
    self.frame.atoms.set_prop("id", Aml.IntPerAtomProperty.from_array(ids))
    self.frame.atoms.set_prop("type", Aml.IntPerAtomProperty.from_array(types))

    IO.dump_frame(frame, simulation.root.joinpath("nonmean").joinpath(f"{step}.trj"))
