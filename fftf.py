#!/usr/bin/env python3

import sys
from pathlib import Path
import openmm as mm
from openmm import unit
from simtk.openmm import app

def save_force(force, filename):
    with open(filename, "w") as f:
        f.write(mm.XmlSerializer.serialize(force))

if __name__ == "__main__":
    topology = app.Topology()
    topology.setPeriodicBoxVectors(
        [
            [ 1, 0, 0 ],
            [ 0, 1, 0 ],
            [ 0, 0, 1 ],

        ]
    )
    ff_file = Path(sys.argv[1]).expanduser().resolve()
    cutoff_nm = float(sys.argv[2])
    ff = app.ForceField(ff_file)

    system = ff.createSystem(topology, nonbondedMethod=app.PME, nonbondedCutoff=cutoff_nm*unit.nanometers, rigidWater=False, removeCMMotion=True)
    forces = system.getForces()

    for i, f in enumerate(forces):
        save_force(f, ff_file.parent.joinpath(ff_file.stem + str(i) + ".xml"))

