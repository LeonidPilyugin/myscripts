#!/usr/bin/env python3

import sys
import matplotlib.pyplot as plt
import numpy as np

plt.style.use("https://raw.githubusercontent.com/LeonidPilyugin/mpl-style/main/simple.mplstyle")

for arg in sys.argv[1:]:
    x = []
    y = []
    p = []
    with open(arg) as f:
        for line in f.readlines()[1:]:
            xx, yy, pp = map(float, line.split(","))
            if yy != 0.0:
                continue
            x.append(xx)
            y.append(yy)
            p.append(pp)

    plt.plot(x, p, label=arg)

if len(sys.argv) > 2:
    plt.legend()

plt.savefig(f"{sys.argv[1]}.png")

