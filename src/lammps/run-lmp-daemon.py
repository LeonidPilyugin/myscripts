#!/usr/bin/env python3

from pathlib import Path
import sys
import time
from tqdm import tqdm
import toml

root = Path(sys.argv[1])
checkpoint_dir = root.joinpath("checkpoint")
trajectory_dir = root.joinpath("trajectory")
log_dir = root.joinpath("log")
src_dir = root.joinpath("src")

with open(root.joinpath("descriptor.toml")) as f:
    data = toml.load(f)["simulation"]

stages = data["stages"]
stage_names = list(stages.keys())

for stage in stages:
    print(f"Stage {stage_names.index(stage) + 1}/{len(stages)}: {stage}")

    t = tqdm(total=stages[stage]["run"])

    max_file = 0
    _max_file = 0

    while not "last.checkpoint" in map(lambda x: x.name, checkpoint_dir.joinpath(stage).iterdir()):
        time.sleep(1)
        if not any(trajectory_dir.joinpath(stage).joinpath(stages[stage]["dumps"][0]["name"]).iterdir()):
            continue
        _max_file = max(map(lambda x: int(x.name.split(".")[0]), trajectory_dir.joinpath(stage).joinpath(stages[stage]["dumps"][0]["name"]).iterdir()))
        _min_file = min(map(lambda x: int(x.name.split(".")[0]), trajectory_dir.joinpath(stage).joinpath(stages[stage]["dumps"][0]["name"]).iterdir()))
        _max_file -= _min_file

        t.update(_max_file - max_file)
        max_file = _max_file

    t.close()
