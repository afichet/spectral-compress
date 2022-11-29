#!/usr/bin/env python3

import os

import matplotlib.pyplot as plt
import numpy as np

path_out  = 'output'

b_12_curves = []
labels = []

for f in os.listdir(path_out):
    f_path = os.path.join(path_out, f)

    with open(f_path) as txt:
        l1 = txt.readline()
        l2 = txt.readline()

        q = [float(v) for v in l1.split()]

        if f.find('_b_12') != -1:
            b_12_curves.append(q)
            labels.append(f)


b_12_curves = np.array(b_12_curves)

b_12_curve_avg = np.average(b_12_curves, axis=0)
x = np.arange(b_12_curves.shape[1])

plt.plot()
# plt.plot(np.arange(b_12_curve_avg.shape[0])[1:], b_12_curve_avg[1:])
# plt.plot(x[1:], b_12_curve_avg[1:])

for q, l in zip(b_12_curves, labels):
    plt.plot(x, q, label=l)

plt.legend()
plt.tight_layout()
plt.show()
