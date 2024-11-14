from gi.repository import Aml
from openmm import unit
import sys

def main(step, simulation, data):
    system = simulation.context.getSystem()
    state = simulation.get_state()
    n = system.getNumParticles()
    positions = state.getPositions(asNumpy=True)
    velocities = state.getVelocities(asNumpy=True)

    if "move_top" in data:
        for i in range(n):
            if positions[i][2] > data["move_top"] * unit.angstrom:
                positions[i][0] += data["move_magnitude"] * unit.angstrom
    else:
        for i in range(n):
            if positions[i][2] < data["move_bottom"] * unit.angstrom:
                positions[i][0] += data["move_magnitude"] * unit.angstrom

    simulation.context.setPositions(positions)
    simulation.context.setVelocities(velocities)

    x = positions[:,0].tolist()
    y = positions[:,1].tolist()
    z = positions[:,2].tolist()
    vx = velocities[:,0].tolist()
    vy = velocities[:,1].tolist()
    vz = velocities[:,2].tolist()
    types = simulation.frame.atoms.get_prop("type").get_arr()
    masses = simulation.frame.atoms.get_prop("mass").get_arr()

    simulation.update_frame(
        n=len(x),
        x=Aml.DoublePerAtomProperty.from_array(x),
        y=Aml.DoublePerAtomProperty.from_array(y),
        z=Aml.DoublePerAtomProperty.from_array(z),
        vx=Aml.DoublePerAtomProperty.from_array(vx),
        vy=Aml.DoublePerAtomProperty.from_array(vy),
        vz=Aml.DoublePerAtomProperty.from_array(vz),
        mass=Aml.DoublePerAtomProperty.from_array(masses),
        type=Aml.IntPerAtomProperty.from_array(types),
    )
