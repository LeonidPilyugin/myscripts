#!/usr/bin/env python3

import re
import argparse
from pathlib import Path

parser = argparse.ArgumentParser("thintraj")
parser.add_argument("-r", type=str, help="Regex to match step number.", default=r"(\d+)")
parser.add_argument("-p", type=str, help="Path to directory with trajectories.", default=".")
parser.add_argument("-k", type=int, help="Keep every files.", default=10)

args = parser.parse_args()

path = Path(args.p)
regex = re.compile(args.r)
keep = args.k

steps = []
size = 0.0
for file in path.glob("*"):
    if not file.is_file():
        continue
    steps.append(int(regex.search(file.name).group(1)))
    size += file.stat().st_size / 2 ** 30

steps.sort()

print(f"Found {len(steps)} files. Total disk usage {size:.2f} G")

keep_steps = steps[::keep]

size = 0.0
for file in path.glob("*"):
    if not file.is_file():
        continue
    if int(regex.search(file.name).group(1)) in keep_steps:
        continue
    size += file.stat().st_size / 2 ** 30
    file.unlink()

print(f"Removed {len(steps) - len(keep_steps)} files. Freed {size:.2f} G")
