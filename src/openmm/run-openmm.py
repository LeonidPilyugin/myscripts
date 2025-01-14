#!/usr/bin/env python3

import sys
import toml
from tqdm import tqdm
from pathlib import Path
import logging
import openmm
import numpy as np
from gi.repository import AmlCore, AmlBasicTypes, AmlParticles, AmlBox, AmlLammpsIo
from addon import OpenmmAddOn, OpenmmSimulation, Logger
from addonloader import AddOnLoader

root = Path(sys.argv[1])

logging.basicConfig(filename=str(root.joinpath("log.log")), level=logging.INFO, format="%(asctime)s: %(message)s")
logger = logging.getLogger()

addon_loader = AddOnLoader(logger=logger)
addon_loader.load()

class Checkpoint:
    root = Path(sys.argv[1])
    chpdir = root.joinpath("checkpoints")
    chpdir.mkdir(exist_ok=True)

    writer_action = AmlLammpsIo.DumpWriter()
    writer_params = AmlLammpsIo.DumpWriterParams()
    writer_params.set_box_id("repr.box")
    writer_params.set_timestep_id("repr.step")
    writer_params.set_particles_id("repr.particles")
    writer_params.set_properties(
        [
            "id",
            "type",
            "mass",
            "x",
            "y",
            "z",
            "vx",
            "vy",
            "vz",
        ]
    )

    reader_action = AmlLammpsIo.DumpReader()
    reader_params = AmlLammpsIo.DumpReaderParams()
    reader_params.set_box_id("box")
    reader_params.set_timestep_id("step")
    reader_params.set_particles_id("particles")
    reader_params.set_properties(
        {
            "id": AmlParticles.Int64PerParticleProperty,
            "type": AmlParticles.Int64PerParticleProperty,
            "mass": AmlParticles.Float64PerParticleProperty,
            "x": AmlParticles.Float64PerParticleProperty,
            "y": AmlParticles.Float64PerParticleProperty,
            "z": AmlParticles.Float64PerParticleProperty,
            "vx": AmlParticles.Float64PerParticleProperty,
            "vy": AmlParticles.Float64PerParticleProperty,
            "vz": AmlParticles.Float64PerParticleProperty,
        }
    )

    @staticmethod
    def get_last_checkpoint_number():
        if list(Checkpoint.chpdir.iterdir()):
            return max([int(file.stem) for file in Checkpoint.chpdir.iterdir()])
        return 0

    @staticmethod
    def load_last_checkpoint(integrator, platform, potentials, potential_data):
        last_checkpoint = Checkpoint.get_last_checkpoint_number()

        if last_checkpoint == 0:
            Checkpoint.reader_params.set_filepath(str(Checkpoint.root.joinpath(f"configuration.lammpsdump")))
        else:
            Checkpoint.reader_params.set_filepath(str(Checkpoint.chpdir.joinpath(f"{last_checkpoint}.lammpsdump")))

        Checkpoint.reader_action.set_params(Checkpoint.reader_params)
        
        # read repr
        repr_obj = AmlCore.DataCollection()

        logging.info(f"Reading lammps dump {Checkpoint.reader_params.get_filepath()}")

        Checkpoint.reader_action.perform(repr_obj)
        repr_obj.set_element("iteration", AmlBasicTypes.Int64.create(last_checkpoint))

        step = 0
        temperature = 0
        kinetic = 0
        potential = 0
        total = 0

        if last_checkpoint > 0:
            logger.info(f"Reading thermo")
            # read thermo
            with open(Checkpoint.chpdir.joinpath(f"{last_checkpoint}.thermo"), "r") as f:
                step, potential, kinetic, total, temperature = map(float, f.read().strip().split())

        repr_obj.set_element("step", AmlBasicTypes.Int64.create(step))
        repr_obj.set_element("thermo.temperature", AmlBasicTypes.Float64.create(temperature))
        repr_obj.set_element("thermo.kinetic", AmlBasicTypes.Float64.create(kinetic))
        repr_obj.set_element("thermo.potential", AmlBasicTypes.Float64.create(potential))
        repr_obj.set_element("thermo.total", AmlBasicTypes.Float64.create(total))

        logger.info("Creating openmm object")

        openmm_obj = OpenmmSimulation()

        particles = repr_obj.get_element("particles")
        
        cell = np.asarray(repr_obj.get_element("box").get_edge().get_arr()).reshape(3, 3) * openmm_obj.unit.length
        positions = np.column_stack((
            np.asarray(particles.get_prop("x").get_arr()),
            np.asarray(particles.get_prop("y").get_arr()),
            np.asarray(particles.get_prop("z").get_arr()),
        )) * openmm_obj.unit.length

        velocities = np.column_stack((
            np.asarray(particles.get_prop("vx").get_arr()),
            np.asarray(particles.get_prop("vy").get_arr()),
            np.asarray(particles.get_prop("vz").get_arr()),
        )) * openmm_obj.unit.velocity

        masses = np.asarray(particles.get_prop("mass").get_arr()) * openmm_obj.unit.mass
        types = np.asarray(particles.get_prop("type").get_arr())

        system = openmm.System()
        system.setDefaultPeriodicBoxVectors(cell[0], cell[1], cell[2])

        for mass in masses:
            system.addParticle(mass)

        for i in range(len(potentials)):
            force = potentials[i]
            data = potential_data[i]

            if hasattr(force, "addParticle"):
                for i in range(particles.get_size()):
                    force.addParticle(*data[str(types[i])])

            system.addForce(force)


        p = openmm.Platform.getPlatformByName(platform["type"])

        context = openmm.Context(system, integrator, p, platform["arguments"])
        context.setPositions(positions)
        context.setVelocities(velocities)

        openmm_obj.context = context
        openmm_obj.integrator = integrator

        return openmm_obj, repr_obj

    @staticmethod
    def dump_checkpoint(data):
        iteration = data.get_element("repr.iteration").get_val()

        logger.info(f"Dumping checkpoint for iteration {iteration}")

        # dump lammps
        Checkpoint.writer_params.set_filepath(str(Checkpoint.chpdir.joinpath(f"{iteration}.lammpsdump")))
        Checkpoint.writer_action.set_params(Checkpoint.writer_params)
        Checkpoint.writer_action.perform(data)

        # dump checkpoint
        openmm_object = data.get_element("openmm")
        with open(Checkpoint.chpdir.joinpath(f"{iteration}.openmmcheckpoint"), "wb") as f:
            f.write(openmm_object.context.createCheckpoint())

        # dump thermo
        step = data.get_element("repr.step").get_val()
        potential = data.get_element("repr.thermo.potential").get_val()
        kinetic = data.get_element("repr.thermo.kinetic").get_val()
        total = data.get_element("repr.thermo.total").get_val()
        temperature = data.get_element("repr.thermo.temperature").get_val()
        with open(Checkpoint.chpdir.joinpath(f"{iteration}.thermo"), "w") as f:
            f.write(f"{step} {potential} {kinetic} {total} {temperature}")


class Load(OpenmmAddOn):
    class Params(AmlCore.ActionParams):
        def __init__(self):
            super().__init__()
            self.integrator = None
            self.platform = None
            self.potentials = None
            self.potential_data = None
            self.sequence = None

        def do_copy(self):
            res = Load.Params()
            res.integrator = self.integrator
            res.platform = self.platform
            res.potentials = self.potentials
            res.potential_data = self.potential_data
            res.sequence = self.sequence
            return res

    class Action(AmlCore.Action):
        def __init__(self):
            super().__init__()

        def do_get_params_error_message(self, params):
            return ""

        def do_perform(self, data):
            params = self.get_params()

            openmm_obj, repr_obj = Checkpoint.load_last_checkpoint(params.integrator, params.platform, params.potentials, params.potential_data)

            logger_obj = Logger()
            logger_obj.logger = logger
            data.set_element("logger", logger_obj)

            data.set_element("openmm", openmm_obj)
            data.set_element("repr", repr_obj)

            data.set_element("filesystem.root", AmlBasicTypes.String.create(sys.argv[1]))
            data.set_element("filesystem.thermofile", AmlBasicTypes.String.create(str(Path(sys.argv[1]).joinpath("thermo.csv"))))

            iteration = repr_obj.get_element("iteration").get_val()

            if iteration == 0:
                logger.info(f"Performing load action sequence")
                for action in params.sequence:
                    logger.info(f"Performing action {action}")
                    action.perform(data)

    def __init__(self):
        self.action = Load.Action()
        self.params = Load.Params()

    def init_dict(self, dictionary):
        # set integrator
        integrator_class = getattr(openmm, dictionary["integrator"]["type"])
        integrator = integrator_class(*dictionary["integrator"]["arguments"])
        self.params.integrator = integrator

        # set platform
        self.params.platform = dictionary["platform"]

        # load potentials
        self.params.potentials = []
        self.params.potential_data = []

        for potential_data in dictionary["potentials"]:
            with open(potential_data["path"], "r") as f:
                self.params.potentials.append(openmm.XmlSerializer.deserialize(f.read()))
                self.params.potential_data.append(potential_data.get("particles", None))

        self.params.sequence = []

        for seq in dictionary.get("sequence", []):
            params = seq.get("params", {})
            classname = seq["type"]

            obj = addon_loader.get_addon_class(classname)()
            obj.init_dict(params)

            self.params.sequence.append(obj.get_action())


        self.action.set_params(self.params)


    def get_action(self):
        return self.action

class Loop(OpenmmAddOn):
    class Params(AmlCore.ActionParams):
        def __init__(self):
            super().__init__()
            self.iterations = None
            self.checkpoint_every = None
            self.sequence = None

        def do_copy(self):
            res = Loop.Params()
            res.iterations = self.iterations
            res.checkpoint_every = self.checkpoint_every
            res.sequence = self.sequence
            return res

    class Action(AmlCore.Action):
        def __init__(self):
            super().__init__()

        def do_get_params_error_message(self, params):
            return ""

        def do_perform(self, data):
            params = self.get_params()

            start = data.get_element("repr.iteration").get_val()
            for i in tqdm(range(start, params.iterations)):
                logger.info(f"Performing iteration {i}")
                data.get_element("repr.iteration").set_val(i)
                for every, action in params.sequence:
                    if i % every == 0:
                        logger.info(f"Performing action {type(action).__name__}")
                        action.perform(data)
                if i % params.checkpoint_every == 0:
                    Checkpoint.dump_checkpoint(data)

    def __init__(self):
        self.action = Loop.Action()
        self.params = Loop.Params()

    def init_dict(self, dictionary):
        self.params.iterations = dictionary["iterations"]
        self.params.checkpoint_every = dictionary["checkpoint_every"]

        self.params.sequence = []

        for seq in dictionary["sequence"]:
            every = seq.get("every", 1)
            params = seq.get("params", {})
            classname = seq["type"]

            obj = addon_loader.get_addon_class(classname)()
            obj.init_dict(params)

            self.params.sequence.append((every, obj.get_action()))

        self.action.set_params(self.params)


    def get_action(self):
        return self.action


if __name__ == "__main__":
    logger.info("Reading descriptor")
    with open(root.joinpath("descriptor.toml")) as f:
        data = toml.load(f)

    load = Load()
    load.init_dict(data["load"])
    load_action = load.get_action()

    loop = Loop()
    loop.init_dict(data["loop"])
    loop_action = loop.get_action()

    data = AmlCore.DataCollection()
    
    load_action.perform(data)
    loop_action.perform(data)

    logger.info("Finished")

