#!/usr/bin/env python3
import sys
import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

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
x, y = np.meshgrid(x, y)
p = np.array(p)
p = p.reshape(x.shape)

surf = ax.plot_surface(x, y, p, linewidth=0, antialiased=False)
plt.show()

