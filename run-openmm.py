#!/usr/bin/env python3

import sys
import toml
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import openmm
import openmm.unit as unit
import os
import numpy as np
import scipy
import scipy.constants
from tqdm import tqdm
import importlib.util

def load_lammps(file):
    with open(file) as f:
        # load number of particles
        for _ in range(3):
            f.readline()
        n_atoms = int(f.readline())
        f.readline()
        
        # load box
        box_coords = []
        for _ in range(3):
            box_coords.append(list(map(float, f.readline().split())))
        
        box = np.array(
            [
                [box_coords[0][1] - box_coords[0][0], 0.0, 0.0],
                [0.0, box_coords[1][1] - box_coords[1][0], 0.0],
                [0.0, 0.0, box_coords[2][1] - box_coords[2][0]],
            ]
        )

        atoms = dict.fromkeys(f.readline().strip("ITEM: ATOMS").split())
        for k in atoms.keys():
            atoms[k] = []

        for line in f.readlines():
            values = line.strip().split()
            for i, val in enumerate(values):
                key = list(atoms.keys())[i]
                val_ = None
                try:
                    val_ = int(val)
                except Exception:
                    val_ = float(val)
                atoms[key].append(val_)

        return atoms, box


def save_lammps(file, box, atoms, step):
    with open(file, "w") as f:
        f.write("ITEM: TIMESTEP\n")
        f.write(f"{step}\n")
        f.write("ITEM: NUMBER OF ATOMS\n")
        f.write(f"{len(atoms['x'])}\n")
        f.write("ITEM: BOX BOUNDS pp pp pp\n")
        f.write(f"0.0 {box[0][0]}\n")
        f.write(f"0.0 {box[1][1]}\n")
        f.write(f"0.0 {box[2][2]}\n")
        f.write("ITEM: ATOMS id ")
        for key in atoms.keys():
            f.write(f"{key} ")
        f.write(f"\n")
        
        for i in range(len(atoms["x"])):
            f.write(f"{i} ")
            for key in atoms.keys():
                f.write(f"{atoms[key][i]} ")
            f.write("\n")

class SimulationData:
    def __init__(self):
        self.forces = []

    def read_ovito(self, filename, length_units=unit.angstroms, time_units=unit.picoseconds, mass_units=unit.atom_mass_units):
        velocity_units = length_units / time_units

        atoms, box = load_lammps(filename)

        self.set_cell(box * length_units)
        self.set_pos(np.column_stack((atoms["x"], atoms["y"], atoms["z"])) * length_units)
        self.set_vel(np.column_stack((atoms["vx"], atoms["vy"], atoms["vz"])) * velocity_units)

        self.masses = atoms["mass"] * mass_units
        self.types = atoms["type"]

    def set_cell(self, cell):
        self.cell = cell
    
    def set_pos(self, pos):
        if hasattr(self, "count"):
            assert self.count == len(pos)

        self.positions = pos
        self.count = len(pos)

    def set_vel(self, vel):
        if hasattr(self, "count"):
            assert self.count == len(vel)

        self.velocities = vel
        self.count = len(vel)

    def set_temp(self, temp):
        self.temperature = temp

    def set_mass(self, mass):
        assert not hasattr(self, "masses")

        self.mass = mass

    def set_integrator(self, integrator):
        assert not hasattr(self, "integrator")

        self.integrator = integrator

    def add_force(self, force):
        self.forces.append(force)

    def add_lj_force(self, sigma, epsilon, cutoff):
        force = openmm.CustomNonbondedForce(
            "4*{epsilon}*(({sigma}/r)^12-({sigma}/r)^6)".format(
                epsilon=epsilon.value_in_unit(unit.kilojoule_per_mole),
                sigma=sigma.value_in_unit(unit.nanometer)
            )
        )
        force.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)
        force.setCutoffDistance(cutoff)
        for _ in range(self.count):
            force.addParticle([])
        self.add_force(force)

    def make_simulation(self, platform, properties):
        try:
            loaded_platform = openmm.Platform.getPlatformByName(platform)
        except openmm.OpenMMException:
            print("Loaded plugins: ", openmm.pluginLoadedLibNames)
            print("Loading errors: ", openmm.Platform.getPluginLoadFailures())
            raise

        system = openmm.System()
        system.setDefaultPeriodicBoxVectors(self.cell[0], self.cell[1], self.cell[2])
        if hasattr(self, "mass"):
            for _ in range(self.count):
                system.addParticle(self.mass)
        else:
            assert hasattr(self, "masses")
            for mass in self.masses:
                system.addParticle(mass)

        for force in self.forces:
            system.addForce(force)
                
        context = openmm.Context(system, self.integrator, loaded_platform, properties)
        
        context.setPositions(self.positions)
        if hasattr(self, 'velocities'):
            context.setVelocities(self.velocities)
        else:
            assert hasattr(self, 'temperature')
            context.setVelocitiesToTemperature(self.temperature)

        simulation = Simulation(context=context, integrator=self.integrator)
        self.set_tainted(True)
        return simulation

    def set_tainted(self, tainted):
        self.tainted = tainted

    def __getattribute__(self, name):
        try:
            if super().__getattribute__('tainted') and name != 'set_tainted':
                raise RuntimeError("Access to a tainted SimulationData")
        except AttributeError:
            pass
        return super().__getattribute__(name)

    def __setattr__(self, name, value) -> None:
        try:
            if super().__getattribute__('tainted') and name != 'tainted':
                raise RuntimeError("Access to a tainted SimulationData")
        except AttributeError:
            pass
        return super().__setattr__(name, value)


class Simulation:
    def __init__(self, context, integrator):
        self.context = context
        self.integrator = integrator
        self.executor = ThreadPoolExecutor(max_workers=1)
        self.running = False

        self.get_state_flags = {'getPositions': True,
                                'getVelocities': True,
                                'enforcePeriodicBox': False,
                                'getForces': True,
                                'getEnergy': True,
                               }
        
        self.masses = np.asarray([context.getSystem().getParticleMass(i).value_in_unit(unit.dalton) / scipy.constants.N_A / 1000 for i in range(context.getSystem().getNumParticles())])

    def get_state(self) -> openmm.State:
        assert not self.running
        return self.context.getState(**self.get_state_flags)

    def step(self, steps):
        assert not self.running
        self.running = True
        self.integrator.step(steps)
        self.running = False
        return self.get_state()

    def step_async(self, steps):
        return self.executor.submit(self.step, (steps))
    
    def dump_ovito(self, state, filename):
        data = ovito.data.DataCollection()

        cell = state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.angstrom)
        positions = state.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
        velocities = state.getVelocities(asNumpy=True).value_in_unit(unit.angstrom / unit.picosecond)

        data_cell = ovito.data.SimulationCell(pbc=(True, True, True))
        data_cell[:, :3] = cell

        data.objects.append(data_cell)

        particles = ovito.data.Particles()
        particles.create_property('Position', data=positions)
        particles.create_property('Velocity', data=velocities)
        particles.create_property('Particle Type', data=types)

        data.objects.append(particles)

        os.makedirs(os.path.dirname(filename), exist_ok=True)

        ovito.io.export_file(
            data,
            filename,
            "lammps/dump",
            columns=[
                "Particle Identifier",
                "Particle Type",
                "Position.X",
                "Position.Y",
                "Position.Z",
                "Velocity.X",
                "Velocity.Y",
                "Velocity.Z"
            ]
        )
        
        
    def mean_next(self, steps, ):
        """Means"""
        
        time_ = self.context.getTime().value_in_unit(unit.second)
        state = self.get_state()
        
        p = state.getPositions(asNumpy=True)
        v = state.getVelocities(asNumpy=True)
        u = t = T = P = 0
        positions = np.zeros_like(p)
        velocities = np.zeros_like(v)
        
        Gs = []

        for _ in range(steps):
            # run 1 step
            state = self.step(1)
            p = state.getPositions(asNumpy=True)
            v = state.getVelocities(asNumpy=True)
            # add parameters to created variables
            positions = np.add(positions, p.value_in_unit(unit.angstrom))
            velocities = np.add(velocities, v.value_in_unit(unit.angstrom / unit.picosecond))
            
            u_ = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
            u += u_
            t_ = state.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)
            t += t_
            
            N_ = len(p)
            
            T_ = (self.masses * np.sum(v.value_in_unit(unit.meter / unit.second) ** 2, axis=1)).sum() / scipy.constants.k / N_ / 3
            T += T_
            
            V_ = state.getPeriodicBoxVolume().value_in_unit(unit.meter ** 3)
            bv_ = state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.meter)
            bv_ = np.asarray([bv_[0][0], bv_[1][1], bv_[2][2]])
            
            p_ = 0 #np.fmod(p.value_in_unit(unit.meter) + (abs(int(np.min(p.value_in_unit(unit.meter))))+1) * bv_, bv_)
            Gs.append((self.masses * np.sum((v.value_in_unit(unit.meter / unit.second) * p_), axis=1)).sum())
        
        time_ = self.context.getTime().value_in_unit(unit.second) - time_
        
        # mean parameters
        positions /= steps
        velocities /= steps
        u /= steps
        t /= steps
        T /= steps
        P = 0 # np.polyfit(np.linspace(0, time_, len(Gs)), Gs, 1)[0] / 3 / V_
        
        return u, t, P, T, positions, velocities, state


    def mean_next_async(self, steps, checkpoint = False):
        return self.executor.submit(self.mean_next, steps, checkpoint)

def dump(therm,
         positions: np.ndarray,
         velocities: np.ndarray,
         u: float,
         t: float,
         P: float,
         T: float,
         step: int,
         cell,
         types,
         atom_file,
    ):
    """Writes dumps of energies and positions"""
    
    therm.write(f"{step},{u},{t},{P},{T}\n")
    therm.flush()

    atoms = {
        "type": types,
        "x": positions[:,0],
        "y": positions[:,1],
        "z": positions[:,2],
    }

    save_lammps(
        atom_file,
        cell,
        atoms,
        step
    )
    
if __name__ == "__main__":
    root = Path(sys.argv[1])
    checkpoint_dir = root.joinpath("checkpoint")
    trajectory_dir = root.joinpath("trajectory")
    log_dir = root.joinpath("log")
    src_dir = root.joinpath("src")

    with open(root.joinpath("descriptor.toml")) as f:
        data = toml.load(f)["simulation"]
    
    # init
    simulation_data = SimulationData()

    # load configuration
    simulation_data .read_ovito(str(src_dir.joinpath("configuration.atom")))

    # load integrator
    integrator_class = getattr(openmm, data["integrator"]["type"])
    integrator = integrator_class(*data["integrator"]["arguments"])
    simulation_data.set_integrator(integrator)

    # load forces
    for force_data in data["potentials"]:
        with open(force_data["path"], "r") as file_force:
            force = openmm.XmlSerializer.deserialize(file_force.read())
        # add particles
        for i in range(simulation_data.count):
            if hasattr(force, "addParticle"):
                force.addParticle(*force_data["particles"][str(simulation_data.types[i])])
        simulation_data.add_force(force)

    types = simulation_data.types

    simulation = simulation_data.make_simulation(data["platform"]["name"], data["platform"]["properties"])

    # compute forces on zero step
    forces = simulation.get_state().getForces(asNumpy=True).value_in_unit(unit.ev / unit.angstrom / unit.mole) / 6.02214076e23
    with open(root.joinpath("forces.csv"), "w") as f:
        f.write("id,fx,fy,fz\n")
        for i, force in enumerate(forces):
            f.write(f"{i},{force[0]},{force[1]},{force[2]}\n")

    
    step = 0
    # load last checkpoint
    if list(checkpoint_dir.iterdir()):
        step = max([int(file.stem) for file in checkpoint_dir.iterdir()])
        checkpoint = str(checkpoint_dir.joinpath(f"{step}.chp"))
        with open(checkpoint, "b+r") as f:
            simulation.context.loadCheckpoint(f.read())

    # create thermo file
    thermo_file = root.joinpath("thermo.csv")
    if step == 0:
        with open(thermo_file, "w") as f:
            f.write("step,u,t,P,T\n")

    iter_steps = data["average_steps"] + data["skip_steps"]
    saved_checkpoints = 0

    # load add-ons
    add_ons = []
    for file in sorted(list(src_dir.joinpath("scripts").iterdir()), key=lambda x: int(x.name.split(".")[0])):
        spec = importlib.util.spec_from_file_location(file.name, file)
        foo = importlib.util.module_from_spec(spec)
        sys.modules[file.name] = foo
        spec.loader.exec_module(foo)
        add_ons.append(foo.main)

    # simulate
    for i in tqdm(range(step, data["steps"], iter_steps)):
        result = simulation.mean_next(data["average_steps"])
        u, t, P, T, p, v, s = result
        with open(thermo_file, "a") as f:
            dump(f, p, v, u, t, P, T, i, s.getPeriodicBoxVectors(asNumpy=True).value_in_unit(openmm.unit.angstrom), types, str(trajectory_dir.joinpath(f"{i}.trj")))
        if data["skip_steps"] > 0:
            simulation.step(data["skip_steps"])

        # save checkpoint
        if data["checkpoint_steps"] > 0 and (i + iter_steps) // data["checkpoint_steps"] >= saved_checkpoints:
            with open(checkpoint_dir.joinpath(f"{i + iter_steps}.chp"), "wb") as ff:
                ff.write(simulation.context.createCheckpoint())
            saved_checkpoints += 1
        for add_on in add_ons:
            simulation, result = add_on(i + iter_steps, simulation, result, data)

    with open(checkpoint_dir.joinpath("last.chp"), "wb") as ff:
        ff.write(simulation.context.createCheckpoint())
