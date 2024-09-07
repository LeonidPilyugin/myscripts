import numpy as np
import openmm
from openmm.openmm import LocalEnergyMinimizer
from openmm import unit
from gi.repository import Aml
import sys

class MinimizationReporter(openmm.MinimizationReporter):
    def report(self, *args):
        #print(args)
        return False

last_inserted = 0

def main(step, simulation, data):
    global last_inserted
    if step != 0 and step - last_inserted < data["insert_xe_every"]:
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

    print(f"Step {step}")

    print("New atom position:", positions[-1])
    print("New atom velocity:", velocities[-1])
    p = [positions[-1] - positions[i] for i in range(len(positions)) if all([abs(positions[-1][j] - positions[i][j]) < 0.5 * unit.nanometer for j in range(3)])]
    print("New atom distances:", p)
    sys.stdout.flush()
    print("After minimization")

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

    # relax
    args = data["emin_args"]
    if "reporter" in args:
        args["reporter"] = MinimizationReporter()
    LocalEnergyMinimizer.minimize(simulation.context, **args)

    # set masses
    for i, m in enumerate(masses):
        system.setParticleMass(i, m)
    last_inserted = step

    state = simulation.get_state()
    positions2 = state.getPositions(asNumpy=True)
    velocities2 = state.getVelocities(asNumpy=True)

    print("New atom position:", positions2[-1])
    print("New atom velocity:", velocities2[-1])
    p = [positions2[-1] - positions2[i] for i in range(len(positions2)) if all([abs(positions2[-1][j] - positions2[i][j]) < 0.5 * unit.nanometer for j in range(3)])]
    print("New atom distances:", p)
    print("\n")
    sys.stdout.flush()

    simulation.context.reinitialize()
    simulation.context.setPositions(positions2)
    simulation.context.setVelocities(velocities2)

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
