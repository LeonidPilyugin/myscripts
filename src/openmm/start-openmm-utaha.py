#!/usr/bin/env python3

import os
import sys
import toml
import json
import shutil
from pathlib import Path

HOME = Path(os.environ.get("SIMULATIONS_HOME", Path.home().joinpath("simulations"))).expanduser()

if __name__ == "__main__":
    assert len(sys.argv) == 2

    with open(sys.argv[1]) as f:
        data = toml.load(f)

    assert "init" in data
    init_data = data.pop("init")

    simulation_id = init_data["id"]
    simulation_path = HOME.joinpath(*simulation_id.split("."))

    assert not simulation_path.exists()
    simulation_path.mkdir(parents=True)
    shutil.copy(init_data["configuration"], simulation_path.joinpath("configuration.lammpsdump"))

    with open(simulation_path.joinpath("descriptor.toml"), "w") as f:
        toml.dump(data, f)

    conda_env = init_data["environment"]

    # start utaha
    descriptor = {
        "taskdata": {
            "alias": f"simulation.{data['id']}",
            "comment": data["comment"],
        },
        "job": {
            "type": "UtahaJobsShellJob",
            "command": [
                "zsh",
                "-c",
                f"source ~/.zshrc && conda activate {conda_env} &&" +
                f"python3 {Path(__file__).parent.joinpath('run-openmm.py')} {simulation_path}"
            ]
        },
        "backend": {
            "type": "UtahaPosixBackendBackend",
        }
    }

    filename = f"/tmp/{os.urandom(32).hex()}.json"

    try:
        with open(filename, "w") as f:
            json.dump(descriptor, f)
        # assert os.system(f"zsh -c 'source ~/.zshrc && utaha --load {filename} && utaha --start simulation.{data['id']}'") == 0
        assert 0 == os.system(f"utaha --load {filename}")
        assert 0 == os.system(f"utaha --start --alias 'simulation.{data['id']}'")
    finally:
        os.unlink(filename)
