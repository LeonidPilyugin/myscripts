#!/usr/bin/env python3

import sys
import toml
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import openmm
import openmm.unit as unit
import ovito.io
import ovito.data
import os
import numpy as np
import scipy
import scipy.constants
from tqdm import tqdm


class SimulationData:
    def __init__(self):
        self.forces = []

    def read_ovito(self, filename, length_units=unit.angstroms, time_units=unit.picoseconds, mass_units=unit.atom_mass_units):
        velocity_units = length_units / time_units

        data = ovito.io.import_file(filename, sort_particles=True).compute()

        self.set_cell(data.cell[:, :3] * length_units)
        self.set_pos(data.particles.positions[...] * length_units)
        self.set_vel(data.particles.velocities[...] * velocity_units)

        self.masses = data.particles.masses[...] * mass_units
        self.types = data.particles.particle_types[...]

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
                print(mass)
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
            
            T_ = (self.masses * np.sum(v.value_in_unit(unit.meter / unit.second) ** 2, axis=1)).sum() / scipy.constants.k / N_ * 2 / 3
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
    ):
    """Writes dumps of energies and positions"""
    
    therm.write(f"{step},{u},{t},{P},{T}\n")
    therm.flush()
    
    data_collection = ovito.data.DataCollection()
    
    # get cell
    data_cell = ovito.data.SimulationCell(pbc=(True, True, True))
    data_cell[:, :3] = cell

    # set cell
    data_collection.objects.append(data_cell)

    # set positions, velocities and types
    particles = ovito.data.Particles()
    particles.create_property("Position", data=positions)
    particles.create_property("Velocity", data=velocities)
    particles.create_property("Particle Type", data=types)
    
    # add data to data object
    data_collection.objects.append(particles)

    # export
    ovito.io.export_file(
        data_collection,
        data["trajectory_template"].format(i=step),
        "lammps/dump",
        columns=[
            "Particle Identifier",
            "Particle Type",
            "Position.X",
            "Position.Y",
            "Position.Z",
            "Velocity.X",
            "Velocity.Y",
            "Velocity.Z",
        ]
    )
    
    return data_collection
    

if __name__ == "__main__":
    root = Path(sys.argv[1])
    checkpoint_dir = root.joinpath("checkpoints")
    trajectory_dir = root.joinpath("trajectory")
    log_dir = root.joinpath("logs")
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
    
    step = 0
    # load last checkpoint
    if checkpoint_dir.iterdir():
        step = max([int(file.name) for file in checkpoint_dir.iterdir()])
        checkpoint = str(checkpoint_dir.joinpath(f"{step}.chp"))
        with open(checkpoint, "b+r") as f:
            simulation.context.loadCheckpoint(f.read())
    
    # create thermo file
    thermo_file = root.joinpath("thermo.csv")
    if step == 0:
        with open(thermo_file, "w") as f:
            f.write("step,u,t,P,T")

    iter_steps = data["average_steps"] + data["skip_steps"]
    saved_checkpoints = 0

    # simulate
    for i in tqdm(range(step, data["steps"] + step, iter_steps)):
        result = simulation.mean_next(data["average_steps"])
        u, t, P, T, p, v, s = result
        with open(thermo_file, "a") as f:
            dump(f, p, v, u, t, P, T, i, s.getPeriodicBoxVectors(asNumpy=True).value_in_unit(openmm.unit.angstrom), types)
        if data["skip_steps"] > 0:
            simulation.step(data["skip_steps"])

        # save checkpoint
        if data["checkpoint_steps"] > 0 and (i + iter_steps) // data["checkpoint_steps"] >= saved_checkpoints:
            with open(checkpoint_dir.joinpath(f"{i + iter_steps}.chp"), "wb") as ff:
                ff.write(simulation.context.createCheckpoint())
            saved_checkpoints += 1

    with open(checkpoint_dir.joinpath("last.chp"), "wb") as ff:
        ff.write(simulation.context.createCheckpoint())
