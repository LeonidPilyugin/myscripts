#!/usr/bin/env python3

import os
import sys
import toml
import shutil
from pathlib import Path

HOME = Path(os.environ.get("SIMULATIONS_HOME", Path.home().joinpath("simulations"))).expanduser()

if __name__ == "__main__":
    with open(sys.argv[1]) as f:
        data = toml.load(f)["simulation"]
    simulation_id = data["id"]
    simulation_path = HOME.joinpath(*simulation_id.split("."))
    input_file = data["input_file"]

    # create directory
    assert not simulation_path.exists()
    simulation_path.mkdir(parents=True)
    simulation_path.joinpath("src").mkdir()
    simulation_path.joinpath("log").mkdir()
    simulation_path.joinpath("trajectory").mkdir()
    simulation_path.joinpath("checkpoint").mkdir()
    shutil.copy(input_file, simulation_path.joinpath("src").joinpath("script.lmp"))
    shutil.copy(sys.argv[1], simulation_path.joinpath("descriptor.toml"))

    conda_env = data["environment"]

    # start spmi
    settings = {
        "task": {
            "id": data["id"],
            "comment": data["comment"],
            "backend": data["spmi"]["backend"],
            "wrapper": {
                "type": "default",
                "command": f"zsh -c 'source ~/.zshrc && conda activate {conda_env} && /usr/bin/env python3 run-lmp.py {simulation_path} && conda deactivate'",
                "mixed_stdout": True,
            },
        }
    }

    filename = f"/tmp/{os.urandom(32).hex()}.toml"

    try:
        with open(filename, "w") as f:
            toml.dump(settings, f)

        assert os.system(f"spmi load {filename}") == 0
        assert os.system(f"spmi start {data['id']}") == 0
    finally:
        os.unlink(filename)

