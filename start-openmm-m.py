#!/usr/bin/env python3

import os
import sys
import toml
import shutil
import subprocess
from pathlib import Path

HOME = Path(os.environ.get("SIMULATIONS_HOME", Path.home().joinpath("simulations"))).expanduser()

if __name__ == "__main__":
    with open(sys.argv[1]) as f:
        data = toml.load(f)["simulation"]

    if not "DeviceIndex" in data["platform"]["properties"]:
        data["platform"]["properties"]["DeviceIndex"] = subprocess.getoutput("select-gpu.py").strip()
    simulation_id = data["id"]
    simulation_path = HOME.joinpath(*simulation_id.split("."))

    # create directory
    assert not simulation_path.exists()
    simulation_path.mkdir(parents=True)
    simulation_path.joinpath("src").mkdir()
    simulation_path.joinpath("scripts").mkdir()
    simulation_path.joinpath("log").mkdir()
    simulation_path.joinpath("trajectory").mkdir()
    simulation_path.joinpath("checkpoint").mkdir()
    shutil.copy(sys.argv[1], simulation_path.joinpath("descriptor.toml"))
    shutil.copy(data["configuration"], simulation_path.joinpath("src").joinpath("configuration.atom"))

    for i in range(len(data["add_ons"])):
        shutil.copy(data["add_ons"][i], simulation_path.joinpath("scripts").joinpath(f"{i}." + Path(data["add_ons"][i]).name))

    conda_env = data["environment"]

    # start spmi
    settings = {
        "task": {
            "id": f"simulation.{data['id']}",
            "comment": data["comment"],
            "backend": data["spmi"]["backend"],
            "wrapper": {
                "type": "default",
                "command": f"zsh -c 'source ~/.zshrc && conda activate {conda_env} && /usr/bin/env python3 {Path(__file__).parent.joinpath('run-openmm-m.py')} {simulation_path} && conda deactivate'",
                "mixed_stdout": True,
            },
        }
    }

    filename = f"/tmp/{os.urandom(32).hex()}.toml"

    try:
        with open(filename, "w") as f:
            toml.dump(settings, f)

        assert os.system(f"zsh -c 'source ~/.zshrc && spmi load {filename} && spmi start simulation.{data['id']}'") == 0
    finally:
        os.unlink(filename)

