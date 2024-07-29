#!/usr/bin/env python3

import sys
import matplotlib.pyplot as plt
import numpy as np

plt.style.use("https://raw.githubusercontent.com/LeonidPilyugin/mpl-style/main/simple.mplstyle")

plt.xlabel("T, K")
plt.ylabel("a, Å")

T = []
DT = []
P = []
DP = []
A = []

for arg in sys.argv[1:]:
    with open(arg) as f:
        for line in f.readlines()[1:]:
            t, dt, p, dp, a, _ = *map(float, line.split(",")[:-1]), ""
            T.append(t)
            DT.append(dt)
            P.append(p)
            DP.append(dp)
            A.append(a)

    plt.plot(T, A, marker=".", label=arg)

if len(sys.argv) > 2:
    plt.legend()

plt.savefig(f"{sys.argv[1]}.png")

