import numpy as np
import openmm
from openmm.openmm import LocalEnergyMinimizer
from openmm import unit
from gi.repository import Aml

last_inserted = 0

def main(step, simulation, data):
    global last_inserted
    if step - last_inserted < data["insert_xe_every"]:
        return

    state = simulation.get_state()
    system = simulation.context.getSystem()

    n = system.getNumParticles()
    positions = state.getPositions(asNumpy=True)
    velocities = state.getVelocities(asNumpy=True)

    # reset masses
    masses = []
    xenons = []
    for i in range(n):
        mass = system.getParticleMass(i)
        if 131 < mass.value_in_unit(unit.amu) < 132:
            xenons.append(i)
        else:
            system.setParticleMass(i, 0.0)
        masses.append(mass)

    # compute COM
    com = np.array([0.0, 0.0, 0.0]) * unit.nanometer
    vel = np.array([0.0, 0.0, 0.0]) * unit.nanometer / unit.picosecond
    for i in xenons:
        com += positions[i]
    com /= len(xenons)

    positions = np.vstack([positions, com]) * unit.nanometer
    velocities = np.vstack([velocities, vel]) * unit.nanometer / unit.picosecond

    # insert Xe
    system.addParticle(masses[xenons[0]])
    masses.append(masses[xenons[0]])

    for i in range(len(system.getForces())):
        force = system.getForce(i)
        if hasattr(force, "addParticle"):
            force.addParticle(*data["potentials"][i]["particles"]["3"])

    simulation.context.reinitialize()
    simulation.context.setPositions(positions)
    simulation.context.setVelocities(velocities)
    simulation.integrator.step(1)
    simulation.context.setPositions(positions)
    simulation.context.setVelocities(velocities)

    # relax
    LocalEnergyMinimizer.minimize(simulation.context, data["emin_tolerance"], data["emin_max_iter"])

    # set masses
    for i in range(n):
        system.setParticleMass(i, masses[i])
    last_inserted = step

    state = simulation.get_state()
    positions2 = state.getPositions(asNumpy=True)
    velocities2 = state.getVelocities(asNumpy=True)

    types = simulation.frame.atoms.get_prop("type").get_arr()
    types.append(3)

    x = positions2[:,0].tolist()
    y = positions2[:,1].tolist()
    z = positions2[:,2].tolist()
    vx = velocities2[:,0].tolist()
    vy = velocities2[:,1].tolist()
    vz = velocities2[:,2].tolist()

    simulation.update_frame(
        n=len(x),
        x=Aml.DoublePerAtomProperty.from_array(x),
        y=Aml.DoublePerAtomProperty.from_array(y),
        z=Aml.DoublePerAtomProperty.from_array(z),
        vx=Aml.DoublePerAtomProperty.from_array(vx),
        vy=Aml.DoublePerAtomProperty.from_array(vy),
        vz=Aml.DoublePerAtomProperty.from_array(vz),
        mass=Aml.DoublePerAtomProperty.from_array([m.value_in_unit(unit.atom_mass_units) for m in masses]),
        type=Aml.IntPerAtomProperty.from_array(types),
    )
