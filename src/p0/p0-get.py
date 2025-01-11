#!/usr/bin/env python3

import sys

def get_isobar(file):
    T = []
    DT = []
    P = []
    DP = []
    P0 = []
    with open(file) as f:
        for line in f.readlines()[1:]:
            t, dt, p, dp, p0, _ = *map(float, line.split(",")[:-1]), ""
            T.append(t)
            DT.append(dt)
            P.append(p)
            DP.append(dp)
            P0.append(p0)

    return T, P0

def predict_a(T, A, temperature):
    indices = []
    if temperature < T[0]:
        indices.extend([0, 1])
    elif temperature >= T[-1]:
        indices.extend([len(T) - 2, len(T) - 1])
    else:
        c = 0
        for t in T:
            if t > temperature:
                break
            c += 1
        indices.extend([c - 1, c])

    dt = T[indices[1]] - T[indices[0]]
    da = A[indices[1]] - A[indices[0]]

    return A[indices[0]] + (temperature - T[indices[0]]) * da / dt


if __name__ == "__main__":
    file = sys.argv[1]
    temperature = float(sys.argv[2])
    T, A = get_isobar(file)
    print(predict_a(T, A, temperature))

