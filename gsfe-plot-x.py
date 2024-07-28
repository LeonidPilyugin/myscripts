#!/usr/bin/env python3

import sys
import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots()

x = []
y = []
p = []
with open(sys.argv[1]) as f:
    for line in f.readlines()[1:]:
        xx, yy, pp = map(float, line.split(","))
        if yy != 0.0:
            continue
        x.append(xx)
        y.append(yy)
        p.append(pp)

plt.plot(x, p)
plt.show()

