#!/usr/bin/env python3

import os
import sys
import toml
import shutil
import subprocess
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

    # start spmi
    settings = {
        "task": {
            "id": f"{init_data['id']}",
            "comment": init_data["comment"],
            "backend": init_data["spmi"]["backend"],
            "wrapper": {
                "type": "default",
                "command": f"zsh -c 'source ~/.zshrc && conda activate {conda_env} && /usr/bin/env python3 {Path(__file__).parent.joinpath('run-openmm.py')} {simulation_path}'",
                "mixed_stdout": True,
            },
        }
    }

    filename = str(simulation_path.joinpath("spmi.toml"))

    with open(filename, "w") as f:
            toml.dump(settings, f)

    os.system(f"zsh -c 'source ~/.zshrc && spmi load {filename} && spmi start {init_data['id']}'")

