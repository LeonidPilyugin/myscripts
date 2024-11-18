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
    positions = state.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
    velocities = state.getVelocities(asNumpy=True).value_in_unit(unit.angstrom / unit.picosecond)

    x = positions[:,0].tolist() * 10
    y = positions[:,1].tolist() * 10
    z = positions[:,2].tolist() * 10
    vx = velocities[:,0].tolist()
    vy = velocities[:,1].tolist()
    vz = velocities[:,2].tolist()
    types = frame.atoms.get_prop("type").get_arr()
    masses = frame.atoms.get_prop("mass").get_arr()
    ids = list(range(1, n+1))

    frame.atoms.set_prop("x", Aml.DoublePerAtomProperty.from_array(x))
    frame.atoms.set_prop("y", Aml.DoublePerAtomProperty.from_array(y))
    frame.atoms.set_prop("z", Aml.DoublePerAtomProperty.from_array(z))
    frame.atoms.set_prop("vx", Aml.DoublePerAtomProperty.from_array(vx))
    frame.atoms.set_prop("vy", Aml.DoublePerAtomProperty.from_array(vy))
    frame.atoms.set_prop("vz", Aml.DoublePerAtomProperty.from_array(vz))
    frame.atoms.set_prop("mass", Aml.DoublePerAtomProperty.from_array(masses))
    frame.atoms.set_prop("id", Aml.IntPerAtomProperty.from_array(ids))
    frame.atoms.set_prop("type", Aml.IntPerAtomProperty.from_array(types))

    IO.dump_frame(frame, str(simulation.root.joinpath("nonmean").joinpath(f"{step}.trj")))
