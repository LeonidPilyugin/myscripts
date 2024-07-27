#!/usr/bin/env python3

import subprocess
from shutil import which
from dataclasses import dataclass

@dataclass
class GpuInfo:
    id: int
    temp: float
    pwr: float
    vram: float


def rocm():
    smi_output = subprocess.check_output("rocm-smi", shell=True).decode("ascii").strip()
    gpus = []
    for line in smi_output.split("\n")[3:-2]:
        arr = line.split()
        gpus.append(
            GpuInfo(
                id = int(arr[0]),
                temp = float(arr[1].strip("c")),
                pwr = float(arr[2].strip("W")),
                vram = float(arr[-2].strip("%")) / 100.0,
            )
        )
    return gpus


def cuda():
    smi_output = subprocess.check_output("nvidia-smi", shell=True).decode("ascii").strip().split("\n")
    gpus = []
    
    for i in range(len(smi_output)):
        line = smi_output[i]
        if "NVIDIA" in line and "CUDA" not in line:
            id = int(line.strip("|").split()[0])
            line = smi_output[i+1]
            line = line.strip("|").split()
            temp = float(line[1].strip("C"))
            pwr = float(line[3].strip("W"))
            vram = float(line[7].strip("MiB")) / float(line[9].strip("MiB"))
            gpus.append(GpuInfo(id=id, temp=temp, pwr=pwr, vram=vram))

    return gpus

def select_gpu(gpus):
    gpus.sort(key=lambda x: x.pwr)
    best_gpu = gpus[0]
    for g in gpus:
        if best_gpu.pwr < g.pwr * 1.1 and g.vram < best_gpu.vram:
            best_gpu = g

    return best_gpu

if __name__ == "__main__":
    gpus = []
    if which("nvidia-smi") is not None:
        gpus.extend(cuda())
    elif which("rocm-smi") is not None:
        gpus.extend(rocm())
    else:
        print("Could not find smi tool")
        exit(1)
    
    print(select_gpu(gpus).id)
