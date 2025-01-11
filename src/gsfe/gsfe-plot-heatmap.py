#!/usr/bin/env python3
import sys
import matplotlib.pyplot as plt
import numpy as np

print("gsfe-plot-heatmap.py file.csv <z> <y> <x>")

x = []
y = []
p = []
with open(sys.argv[1]) as f:
    for line in f.readlines()[1:]:
        xx, yy, pp = map(float, line.split(","))
        x.append(xx)
        y.append(yy)
        p.append(pp)

x = list(dict.fromkeys(x))
y = list(dict.fromkeys(y))

x = np.array(x)
y = np.array(y)
p = np.array(p)

p = -p

p = p.reshape((x.size, y.size))

plt.imshow(p, cmap="Spectral", interpolation="nearest")

# remove ticks
plt.xticks([], [])
plt.yticks([], [])

plt.title(f"{{{sys.argv[2]}}}")
plt.xlabel(f"<{sys.argv[3]}>")
plt.ylabel(f"<{sys.argv[4]}>")

plt.savefig(f"{sys.argv[1]}.png")
