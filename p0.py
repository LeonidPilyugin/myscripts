#!/usr/bin/env python3

import os
import sys
import toml
import subprocess
from time import sleep
from random import randint
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
    processes: int = 4

@dataclass
class GDParams:
    iterations: int = 20
    parameters: int = 0
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

    temperatures = { i: [] for i in range(params.processes) }
    pressures = { i: [] for i in range(params.processes) }
    processes = []

    try:
        for i in range(params.processes):
            mscript = script.format(random_seed=randint(1, 2 ** 30))
            with open(filename + str(i), "w") as f:
                f.write(mscript)

            command = f"{params.lammps_executable} -in {filename}{i} -l {filename}{i}.log",
            processes.append(
                subprocess.Popen(
                    command,
                    shell=True,
                    stdout=sys.stderr,
                    stderr=sys.stderr,
                )
            )
            sleep(1)
        for i in range(params.processes):
            processes[i].wait()
            assert processes[i].returncode == 0

            with open(f"{filename}{i}.log") as f:
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
                            temperatures[i].append(temp)
                            pressures[i].append(pres)
                        except Exception:
                            continue
    finally:
        os.system("rm {filename}*")
        pass

    ts = [ np.mean(x) for x in temperatures.values() ]
    dts = [ np.std(x) for x in temperatures.values() ]
    ps = [ np.mean(x) for x in pressures.values() ]
    dps = [ np.std(x) for x in pressures.values() ]

    best = ps.index(min([abs(p) for p in ps]))

    temperature, current_pressure = (ts[best], dts[best]), (ps[best], dps[best])

    print(f"{params.parameters}: T = {temperature[0]} +/- {temperature[1]}, P = {current_pressure[0]} +/- {current_pressure[1]}")
    sys.stdout.flush()

    return temperature, current_pressure

def optimize(params: Params, current_params=None):

    if current_params is None:
        current_params = np.array([5.0] * parameters.gd_params.parameters)

    # create start parameters
    params_history = [current_params, current_params * 0.9]
    pressure_history = []
    temperature_history = []

    def lammps(p):
        params.lammps_params.parameters = p
        t, p = lammps_run(params.lammps_params)
        pressure_history.append(p)
        temperature_history.append(t)

    for p in params_history:
        lammps(p)

    # get best params
    while len(pressure_history) < 3 or all(
        [
            abs(pressure_history[-1][0]) < abs(pressure_history[-2][0]),
            abs(pressure_history[-1][0]) > pressure_history[-1][1],
        ]
    ):
        params_history.append(
            params_history[-1] - pressure_history[-1][0] * (params_history[-1] - params_history[-2]) / (pressure_history[-1][0] - pressure_history[-2][0]) * params.gd_params.learning_rate
        )

        lammps(params_history[-1])

    mp = [x[0] for x in pressure_history]
    best = mp.index(min(mp))

    return params_history[best], *pressure_history[best], *temperature_history[best]



def process(params: Params):
    param_list = []
    parameters = None
    for temperature in np.arange(params.t_start, params.t_stop, params.t_step):
        print(f"Processing {temperature}")
        params.lammps_params.temperature = temperature
        parameters, pres, dpres, temp, dtemp = optimize(params, parameters)
        param_list.append((parameters, temp, dtemp, pres, dpres))

    return param_list

def print_parameters(param_list: list, f):
    first = param_list[0]
    f.write("T,DT,P,DP,", end="")
    for i in range(len(first[0])):
        f.write(f"P{i},", end="")
    f.write("\n")
    for param in param_list:
        f.write(f"{param[1]},{param[2]},{param[3]},{param[4]},", end="")
        for p in param[0]:
            f.write(f"{p},")
        f.write("\n")


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
    with open(sys.argv[2], "w") as f:
        print_parameters(param_list, f)
  
