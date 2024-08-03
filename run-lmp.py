#!/usr/bin/env python3

import sys
import toml
import subprocess
from pathlib import Path


IGNORE = [
    "dump",
    "undump",
    "fix",
    "unfix",
    "compute",
    "uncompute",
    "thermo",
    "thermo_style",
    "thermo_modify",
    "pair_style",
    "pair_coeff",
]

root = Path(sys.argv[1])
checkpoint_dir = root.joinpath("checkpoint")
trajectory_dir = root.joinpath("trajectory")
log_dir = root.joinpath("log")
src_dir = root.joinpath("src")

with open(root.joinpath("descriptor.toml")) as f:
    data = toml.load(f)["simulation"]

data["macros"]["__ROOT"] = str(root)

stages = data["stages"]
stage_names = list(stages.keys())

for name in stage_names:
    p = checkpoint_dir.joinpath(name)
    if not p.exists(): p.mkdir()
    p = trajectory_dir.joinpath(name)
    if not p.exists(): p.mkdir()

def find_checkpoint():
    stage = None
    checkpoint = None
    steps = None

    for name in stage_names:
        directory = checkpoint_dir.joinpath(name)
        if "last.checkpoint" in map(lambda x: x.name, directory.iterdir()):
            stage = name
            checkpoint = list(filter(lambda x: x.name == "last.checkpoint", directory.iterdir()))[0]
            steps = 0
            continue
        if not any(directory.iterdir()):
            break

        stage = name
        max_num = max(map(lambda x: int(x.name.split('.')[0]), list(directory.iterdir())))
        checkpoint = directory.joinpath(f"{max_num}.checkpoint")
        steps = stages[stage]["run"] - max_num
        
        for stg in stage_names[:stage_names.index(name)]:
            steps += stages[stg]["run"]

        for trj_file in trajectory_dir.joinpath(name).iterdir():
            if int(trj_file.name.split(".")[0]) > max_num:
                trj_file.unlink()
        break

    return (stage, checkpoint, steps)

def modify_script(stage, checkpoint, steps):
    script = src_dir.joinpath("script.lmp")
    modified_script = src_dir.joinpath("modified-script.lmp")

    stage_n = 0

    with open(script) as input_file, open(modified_script, "w") as output_file:
        if not stage is None:
            output_file.write(f"read_restart '{checkpoint}'\n")

        for line in input_file.readlines():
            if stage_n >= len(stages):
                output_file.write(line)
                continue
            stage_name = stage_names[stage_n]
            if not line.strip().split() or line.strip().split()[0] != "run":
                line = line.replace("__MEAN", str(stages[stage_name]["mean"]))
                line = line.replace("__DUMP_EVERY", str(stages[stage_name]["dump"]))
                line = line.replace("__DUMP_PATH", str(trajectory_dir.joinpath(stage_name)) + "/*.trj")
                line = line.replace("__TIMESTEP", str(stages[stage_name]["timestep"]))
                line = line.replace("__THERMO_EVERY", str(stages[stage_name]["thermo"]))
                line = line.replace("__THERMO_MEAN", str(stages[stage_name]["thermo_mean"]))
                
                if not stage is None and stage_names.index(stage) >= stage_names.index(stage_name):
                    if len(line.strip().split()) == 0 or line.strip().split()[0] not in IGNORE:
                        line = "# " + line

                output_file.write(line)
            else:
                if stage is None or stage_names.index(stage) <= stage_names.index(stage_name):
                    if stages[stage_name]['checkpoint'] > 0:
                        output_file.write(f"restart {stages[stage_name]['checkpoint']} '{checkpoint_dir.joinpath(stage_name)}/*.checkpoint'\n")

                    output_file.write(f"timestep {stages[stage_name]['timestep']}\n")

                    output_file.write(f"thermo {stages[stage_name]["thermo"]}\n")

                    if stage == stage_name:
                        output_file.write(f"run {steps}\n")
                    else:
                        output_file.write(f"run {stages[stage_name]['run']}\n")

                    output_file.write(f"write_restart '{checkpoint_dir.joinpath(stage_name)}/last.checkpoint'\n")
                stage_n += 1

    with open(modified_script) as read_file:
        contents = read_file.read()

    for key, value in data["macros"].items():
        contents = contents.replace(key, value)

    with open(modified_script, "w") as write_file:
        write_file.write(contents)

    return modified_script


def get_log():
    max_log = 0 if not any(log_dir.iterdir()) else max(map(lambda x: int(x.name.split('.')[0]), log_dir.iterdir())) + 1
    return f"{str(log_dir)}/{max_log}.log"


def start(script):
    pb = subprocess.Popen(
        f"/usr/bin/env python3 {Path(__file__).parent.joinpath('run-lmp-daemon.py')} '{sys.argv[1]}'",
        shell=True,
    )

    vargs = ""
    for key, value in data["variables"].items():
        vargs += f"-var {key} {value} "

    p = subprocess.Popen(
        f"{data['executable']} -nonbuf -l {get_log()} -in '{script}' {vargs} {data['flags']}",
        shell=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    p.wait()
    if p.returncode == 0:
        pb.wait()
    else:
        pb.kill()
    sys.exit(p.returncode)

start(modify_script(*find_checkpoint()))

