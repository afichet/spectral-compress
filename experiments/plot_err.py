#!/usr/bin/env python3

import os

import matplotlib.pyplot as plt
import numpy as np

path_out  = 'output'

b_bits = []
b_err = []

u_bits = []
u_err = []

utb_bits = []
utb_err = []

up_bits = []
up_err = []

for f in os.listdir(path_out):
    f_path = os.path.join(path_out, f)

    with open(f_path) as txt:
        l1 = txt.readline()
        l2 = txt.readline()

        q = np.sum([float(v) for v in l1.split()])
        e = float(l2)
        if e > 1:
            print(f_path)
        else:
            if f.find('_b_') != -1:
                b_bits.append(q)
                b_err.append(e)
            elif f.find('_u_') != -1:
                u_bits.append(q)
                u_err.append(e)
            elif f.find('_utb_') != -1:
                utb_bits.append(q)
                utb_err.append(e)
            elif f.find('_up_') != -1:
                up_bits.append(q)
                up_err.append(e)


plt.plot()
plt.xlabel('Bits per pixel')
plt.ylabel('Error')

pt_sz = 8

plt.scatter(b_bits, b_err, s=pt_sz, label='Bounded')
plt.scatter(u_bits, u_err, s=pt_sz, label='Unbounded')
plt.scatter(utb_bits, utb_err, s=pt_sz, label='Unbounded to bounded')
plt.scatter(up_bits, up_err, s=pt_sz, label='Upperbound')

plt.legend()
plt.tight_layout()
plt.show()
