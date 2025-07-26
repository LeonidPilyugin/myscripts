#!/usr/bin/env python3

import os
import sys
import toml
import json
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
    shutil.copy(Path(input_file).expanduser(), simulation_path.joinpath("src").joinpath("script.lmp"))
    shutil.copy(Path(sys.argv[1]).expanduser(), simulation_path.joinpath("descriptor.toml"))

    conda_env = data["environment"]

    # start utaha

    descriptor = {
        "taskdata": {
            "alias": f"simulation.{data['id']}",
            "comment": data["comment"],
        },
        "wrapper": {
            "type": "UtahaDefaultWrapperWrapper",
            "command": [
                "zsh",
                "-c",
                f"source ~/.zshrc && conda activate {conda_env} &&" +
                f"python3 {Path(__file__).parent.joinpath('run-lmp.py')} {simulation_path}"
            ]
        },
        "backend": {
            "type": "UtahaScreenBackendBackend",
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
