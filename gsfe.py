#!/usr/bin/env python3

import os
import sys
import toml
import subprocess
from dataclasses import dataclass, asdict

LAMMPS_FORMAT = r"""
boundary p p f
units metal
atom_style {atom_style}

{create_box}

{load_potential}

compute PE all pe 

region TOP block $(xlo) $(xhi) $(ylo) $(yhi) $((zlo + zhi) / 2) $(zhi) units box
region BOTTOM block $(xlo) $(xhi) $(ylo) $(yhi) $(zlo) $((zlo + zhi) / 2) units box

label loop_x
variable step_x loop 0 {steps}

label loop_y
variable step_y loop 0 {steps}

variable dx equal $(v_step_x * {x_displ} / {steps})
variable dy equal $(v_step_y * {y_displ} / {steps})

{create_atoms}

group TOP_GROUP region TOP 
group BOTTOM_GROUP region BOTTOM 

displace_atoms TOP_GROUP move $(v_dx) $(v_dy) 0 units box 
fix SETFORCE_FIX all setforce 0.0 0.0 NULL
minimize {minimize_params}

unfix SETFORCE_FIX
print "GSFE_DATA: $(v_dx),$(v_dy),$(c_PE)"
delete_atoms group all

next step_y
jump SELF loop_y
next step_x
jump SELF loop_x
"""
@dataclass
class Params:
    lammps_executable: str = ""
    atom_style: str = "atomic"
    create_box: str = ""
    load_potential: str = ""
    create_atoms: str = ""
    steps: int = 0 
    x_displ: float = 0.0
    y_displ: float = 0.0
    minimize_params: str = ""

@dataclass
class Point:
    p: float = 0.0
    x: float = 0.0
    y: float = 0.0


def lammps_run(params: Params):
    script = LAMMPS_FORMAT.format(
        **asdict(params),
    )

    filename = f"/tmp/{os.urandom(32).hex()}.lmp"

    points = []

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
            for line in f.readlines():
                if line.strip().startswith("GSFE_DATA:"):
                    data = list(map(float, line.strip().split()[1].split(",")))
                    points.append(
                        Point(
                            x = data[0],
                            y = data[1],
                            p = data[2],
                        )
                    )

    finally:
        os.unlink(filename)
        os.unlink(f"{filename}.log")
        pass

    return points


def print_points(points):
    e0 = 0.0
    for point in points:
        if point.x == point.y == 0.0:
            e0 = point.p
            break
    for point in points:
        point.p -= e0

    print("x,y,e")
    for point in points:
        print(f"{point.x},{point.y},{point.p}")

if __name__ == "__main__":
    with open(sys.argv[1]) as f:
        data = toml.load(f)

    parameters = Params(
        lammps_executable = data["lammps_executable"],
        atom_style = data["atom_style"],
        create_box = data["create_box"],
        load_potential = data["load_potential"],
        create_atoms = data["create_atoms"],
        steps = data["steps"],
        x_displ = data["x_displ"],
        y_displ = data["y_displ"],
        minimize_params = data["minimize_params"],
    )

    points = lammps_run(parameters)
    print_points(points)
  
