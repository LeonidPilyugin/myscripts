#!/usr/bin/env python3

import os
import sys
import toml
import subprocess
from dataclasses import dataclass, field, asdict
import numpy as np

LAMMPS_FORMAT = r"""
# ========== begin ==========
boundary p p p
units metal
atom_style {atom_style}

# ========== simulation box ==========
{create_box}

# ========== create reference system ==========
{create_atoms}

# ========== load potential ==========
{load_potential}

# ========== set velocity ==========
velocity all create {temperature} {random_seed}
fix NVT_FIX all nvt temp {temperature} {temperature} {nvt_ng_param}
run {relax_steps}
unfix NVT_FIX

# compute temperature, pressure and potential energy
compute TEMP all temp
compute PRESSURE all pressure TEMP
compute PE all pe
# compute average values
fix TEMP_FIX all ave/time 1 {mean} {frame} c_TEMP
fix PRESSURE_FIX all ave/time 1 {mean} {frame} c_PRESSURE
# set thermo frequency
thermo {frame}
thermo_style custom step f_TEMP_FIX f_PRESSURE_FIX

# ========== start ==========
fix FIX_NVE all nve
timestep {timestep}
run {steps}
"""
@dataclass
class LammpsParams:
    temperature: int = 0
    parameters: np.ndarray = field(default_factory=lambda: np.array([]))
    create_box: str = ""
    load_potential: str = ""
    create_atoms: str = ""
    atom_style: str = "atomic"
    steps: int = 50000
    relax_steps: int = 50000
    mean: int = 1000
    frame: int = 10000
    lammps_executable: str = ""
    random_seed: int = 10
    nvt_ng_param: float = 0.3
    timestep: float = 1e-3

@dataclass
class GDParams:
    iterations: int = 20
    parameters: int = 0
    pressure_error: int = 50
    temperature_error: int = 20
    learning_rate = 0.8

@dataclass
class Params:
    lammps_params: LammpsParams = field(default_factory=lambda: LammpsParams())
    gd_params: GDParams = field(default_factory=lambda: GDParams())
    t_start: float = 0
    t_stop: float = 0
    t_step: float = 0


def lammps_run(params: LammpsParams):
    script = LAMMPS_FORMAT.format(
        **asdict(params),
    )

    script = script.format(
        **{ f"P{i}": params.parameters[i] for i in range(len(params.parameters)) }
    )

    filename = f"/tmp/{os.urandom(32).hex()}.lmp"

    temperatures = []
    pressures = []

    try:
        with open(filename, "w") as f:
            f.write(script)

        command = f"{params.lammps_executable} -in {filename} -l {filename}.log",
        p = subprocess.run(
            command,
            shell=True,
            stdout=sys.stderr,
            stderr=sys.stderr,
        )
        assert p.returncode == 0

        with open(f"{filename}.log") as f:
            read_thermo = False
            for line in f.readlines():
                if line.startswith("Per MPI"):
                    read_thermo = True
                    continue
                if line.startswith("Loop time"):
                    read_thermo = False
                    continue
                if read_thermo:
                    try:
                        _, temp, pres = map(float, line.split())
                        temperatures.append(temp)
                        pressures.append(pres)
                    except Exception:
                        continue
    finally:
        os.unlink(filename)
        os.unlink(f"{filename}.log")
        pass

    return (np.mean(temperatures[1:]), np.std(temperatures[1:])), (np.mean(pressures[1:]), np.std(pressures[1:]))


def gradient_descent(params: Params, current_params=None):
    if current_params is None:
        current_params = np.array([5.0] * parameters.gd_params.parameters)
    previous_params = np.array([0.0] * parameters.gd_params.parameters)
    current_pressure = None
    previous_pressure = None
    temperature = None

    delta = 10

    while current_pressure is None or all(
        [
            abs(current_pressure[0]) >= current_pressure[1],
            delta >= 1e-6
        ]
    ):
        params.lammps_params.parameters = current_params
        temperature, current_pressure = lammps_run(params.lammps_params)

        print(f"{current_params}: T = {temperature[0]} +/- {temperature[1]}, P = {current_pressure[0]} +/- {current_pressure[1]}")
        sys.stdout.flush()

        if previous_pressure is None:
            previous_pressure = current_pressure[0]
            previous_params = current_params
            current_params = current_params * 0.9
            continue


        new_params = np.abs(
            current_params - current_pressure[0] * (current_params - previous_params) / (current_pressure[0] - previous_pressure) * params.gd_params.learning_rate
        )

        delta = np.max(np.abs((new_params - current_params) / current_params))

        previous_pressure = current_pressure[0]
        previous_params = current_params
        current_params = new_params

    return current_params, current_pressure[0], current_pressure[1], temperature[0], temperature[1]



def process(params: Params):
    param_list = []
    parameters = None
    for temperature in np.arange(params.t_start, params.t_stop, params.t_step):
        print(f"Processing {temperature}")
        params.lammps_params.temperature = temperature
        parameters, pres, dpres, temp, dtemp = gradient_descent(params, parameters)
        param_list.append((parameters, temp, dtemp, pres, dpres))

    return param_list

def print_parameters(param_list: list):
    first = param_list[0]
    print("T,DT,P,DP,", end="")
    for i in range(len(first[0])):
        print(f"P{i},", end="")
    print()
    for param in param_list:
        print(f"{param[1]},{param[2]},{param[3]},{param[4]}", end="")
        for p in param[0]:
            print(f"{p},")


if __name__ == "__main__":
    with open(sys.argv[1]) as f:
        data = toml.load(f)

    parameters = Params()
    parameters.t_start = float(data["t_start"])
    parameters.t_stop = float(data["t_stop"])
    parameters.t_step = int(data["t_step"])
    parameters.gd_params.parameters = int(data["parameters"])
    parameters.lammps_params.create_box = data["create_box"]
    parameters.lammps_params.atom_style = data["atom_style"]
    parameters.lammps_params.load_potential = data["load_potential"]
    parameters.lammps_params.create_atoms = data["create_atoms"]
    parameters.lammps_params.lammps_executable = data["lammps_executable"]

    param_list = process(parameters)
    print_parameters(param_list)
  
