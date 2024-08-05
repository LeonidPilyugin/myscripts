#!/usr/bin/env python3

import sys
from pathlib import Path
import numpy as np

KJ_NM_TO_EV_A = 0.0010364269656262175

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
    print(np.max(np.abs(delta)))
    print(np.mean(np.abs(delta)))

