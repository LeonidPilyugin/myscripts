#!/usr/bin/env python3

import sys
import math
from pathlib import Path
import numpy as np

KJ_NM_TO_EV_A = 1.0367955273190636e-05

def load_csv(file):
    res = []
    with open(file) as f:
        for line in f.readlines()[1:]:
            res.append(list(map(float, line.split(",")[1:])))

    return np.array(res) * KJ_NM_TO_EV_A

def load_lammps(file):
    res = []
    with open(file) as f:
        for line in f.readlines()[9:]:
            res.append(list(map(float, line.split()[1:])))

    return np.array(res)



if __name__ == "__main__":
    files = list(map(Path, sys.argv[1:]))
    forces = []
    for file in files:
        if file.suffix == ".csv":
            forces.append(load_csv(file))
        else:
            forces.append(load_lammps(file))

    delta = forces[0] - forces[1]
    delta /= forces[0]
    total = [ math.hypot(*d) for d in delta ]
    print(f"Max reliative error: {max(total)}")
    print(f"Mean reliative error: {np.mean(total)}")

