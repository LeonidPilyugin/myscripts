#!/usr/bin/env python3

import sys
import toml
from pathlib import Path
import openmm
import openmm.unit as unit
import os
import numpy as np
import scipy
import scipy.constants
from tqdm import tqdm
import importlib.util
from gi.repository import Aml

IO = Aml.Io.create(Aml.LammpsTextDumpParser())
LENGTH_UNITS = unit.angstroms
TIME_UNITS = unit.picoseconds
VELOCITY_UNITS = LENGTH_UNITS / TIME_UNITS
MASS_UNITS = unit.atom_mass_units

class SimulationData:
    def __init__(self):
        self.forces = []

    def load(self, filename):
        self.frame = IO.load_frame(filename)

        self.set_cell(np.asarray(self.frame.box.get_edge().get_arr()).reshape((3, 3)) * LENGTH_UNITS)

        self.set_pos(np.column_stack(
            (
                np.asarray(self.frame.atoms.get_prop("x").get_arr()),
                np.asarray(self.frame.atoms.get_prop("y").get_arr()),
                np.asarray(self.frame.atoms.get_prop("z").get_arr()),
            )
        ) * LENGTH_UNITS)

        self.set_vel(np.column_stack(
            (
                np.asarray(self.frame.atoms.get_prop("vx").get_arr()),
                np.asarray(self.frame.atoms.get_prop("vy").get_arr()),
                np.asarray(self.frame.atoms.get_prop("vz").get_arr()),
            )
        ) * VELOCITY_UNITS)

        self.masses = np.asarray(self.frame.atoms.get_prop("mass").get_arr()) * MASS_UNITS
        self.types = np.asarray(self.frame.atoms.get_prop("type").get_arr())

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

        simulation = Simulation(context, self.integrator, self.frame)
        return simulation


class Simulation:
    def __init__(self, context, integrator, frame):
        self.frame = frame
        self.context = context
        self.integrator = integrator
        self.get_state_flags = {'getPositions': True,
                                'getVelocities': True,
                                'enforcePeriodicBox': False,
                                'getForces': True,
                                'getEnergy': True,
                               }

    def get_state(self) -> openmm.State:
        return self.context.getState(**self.get_state_flags)

    def skip_steps(self, steps):
        self.integrator.step(steps)
        
    def mean_next(self, steps):
        """Means"""
        
        time_ = self.context.getTime().value_in_unit(unit.second)
        state = self.get_state()
        
        p = state.getPositions(asNumpy=True)
        v = state.getVelocities(asNumpy=True)
        u = t = T = 0
        positions = np.zeros_like(p)
        velocities = np.zeros_like(v)

        masses = np.asarray(self.frame.atoms.get_prop("mass").get_arr())
        masses /= scipy.constants.N_A * 1000
        
        for _ in range(steps):
            # run 1 step
            self.integrator.step(1)
            state = self.get_state()
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
            
            T_ = (masses * np.sum(v.value_in_unit(unit.meter / unit.second) ** 2, axis=1)).sum() / scipy.constants.k / N_ / 3
            T += T_
            
            bv_ = state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.meter)
            bv_ = np.asarray([bv_[0][0], bv_[1][1], bv_[2][2]])
            
        time_ = self.context.getTime().value_in_unit(unit.second) - time_
        
        # mean parameters
        positions /= steps
        velocities /= steps
        u /= steps
        t /= steps
        T /= steps

        self.frame.set_prop("u", Aml.FrameProperty.create(u))
        self.frame.set_prop("t", Aml.FrameProperty.create(t))
        self.frame.set_prop("T", Aml.FrameProperty.create(T))

        x = positions[:,0].tolist()
        y = positions[:,1].tolist()
        z = positions[:,2].tolist()
        vx = velocities[:,0].tolist()
        vy = velocities[:,1].tolist()
        vz = velocities[:,2].tolist()

        self.update_frame(
            n=len(x),
            x=Aml.DoublePerAtomProperty.from_array(x),
            y=Aml.DoublePerAtomProperty.from_array(y),
            z=Aml.DoublePerAtomProperty.from_array(z),
            vx=Aml.DoublePerAtomProperty.from_array(vx),
            vy=Aml.DoublePerAtomProperty.from_array(vy),
            vz=Aml.DoublePerAtomProperty.from_array(vz),
            mass=self.frame.atoms.get_prop("mass"),
            type=self.frame.atoms.get_prop("type"),
        )

    def update_frame(self, *args, **props):
        self.frame.atoms = Aml.Atoms.sized(props["n"])
        self.frame.atoms.set_prop("x", props["x"])
        self.frame.atoms.set_prop("y", props["y"])
        self.frame.atoms.set_prop("z", props["z"])
        self.frame.atoms.set_prop("vx", props["vx"])
        self.frame.atoms.set_prop("vy", props["vy"])
        self.frame.atoms.set_prop("vz", props["vz"])
        self.frame.atoms.set_prop("mass", props["mass"])
        self.frame.atoms.set_prop("type", props["type"])
        
    def dump(self, step, thermo_file, trj_file):
        u = self.frame.get_prop("u").data
        t = self.frame.get_prop("t").data
        T = self.frame.get_prop("T").data

        with open(thermo_file, "a") as f:
            f.write(f"{step},{u+t},{u},{t},{T}\n")
            f.flush()

        IO.dump_frame(self.frame, trj_file)


if __name__ == "__main__":
    root = Path(sys.argv[1])
    checkpoint_dir = root.joinpath("checkpoint")
    trajectory_dir = root.joinpath("trajectory")
    log_dir = root.joinpath("log")
    src_dir = root.joinpath("src")

    with open(root.joinpath("descriptor.toml")) as f:
        data = toml.load(f)["simulation"]

    simulation_data = SimulationData()

    if list(checkpoint_dir.iterdir()):
        step = max([int(file.stem) for file in checkpoint_dir.iterdir()])
        simulation_data.load(str(trajectory_dir.joinpath(f"{step}.trj")))
    else:
        step = 0
        simulation_data.load(str(src_dir.joinpath("configuration.atom")))

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
            f.write("step,e,u,t,T\n")
    else:
        with open(thermo_file, "r") as f:
            lines = []
            st = 0
            line = f.readline()
            while st != step:
                lines.append(line)
                line = f.readline()
                st = int(line.split()[0])

        with open(thermo_file, "w") as f:
            for line in lines:
                f.write(line)

    iter_steps = data["average_steps"] + data["skip_steps"]
    saved_checkpoints = 0

    # load add-ons
    add_ons = []
    for file in sorted(list(filter(lambda x: x.is_file(), root.joinpath("scripts").iterdir())), key=lambda x: int(x.name.split(".")[0])):
        spec = importlib.util.spec_from_file_location(file.name, file)
        foo = importlib.util.module_from_spec(spec)
        sys.modules[file.name] = foo
        spec.loader.exec_module(foo)
        add_ons.append(foo.main)

    # simulate
    for i in tqdm(range(step, data["steps"], iter_steps)):
        simulation.mean_next(data["average_steps"])
        simulation.dump(i, thermo_file, str(trajectory_dir.joinpath(f"{i}.trj")))

        for add_on in add_ons:
            add_on(i, simulation, data)

        if data["skip_steps"] > 0:
            simulation.skip_steps(data["skip_steps"])

        # save checkpoint
        if data["checkpoint_steps"] > 0 and (i + iter_steps) // data["checkpoint_steps"] >= saved_checkpoints:
            with open(checkpoint_dir.joinpath(f"{i + iter_steps}.chp"), "wb") as ff:
                ff.write(simulation.context.createCheckpoint())
            saved_checkpoints += 1

    with open(checkpoint_dir.joinpath("last.chp"), "wb") as ff:
        ff.write(simulation.context.createCheckpoint())
