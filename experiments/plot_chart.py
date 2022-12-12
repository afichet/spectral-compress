#!/usr/bin/env python3

import matplotlib
import matplotlib.pyplot as plt
import numpy as np


use_latex = True

if use_latex:
    matplotlib.use("pgf")
    matplotlib.rcParams.update({
        "pgf.texsystem": "pdflatex",
        'font.family': 'serif',
        'text.usetex': True,
        'pgf.rcfonts': False,
    })

# Execution parameters
dataset = 'CAVE/chart_and_stuffed_toy_ms'
dataset = 'bonn/VereinsStoffe_Seite_nach_Notenschluessel_rotkrone_aniso_medium_gloss__diffuse'

u_bits = []
u_err = []

ub_bits = []
ub_err = []

techniques = [
    # ('bounded', 'b'),
    ('unbounded', 'u', [], [], [], 'tab:blue'),
    ('unbounded to bounded', 'utb', [], [], [], 'tab:orange'),
    ('upperbound', 'up', [], [], [], 'tab:green')
]

n_bits = [6, 7, 9, 10, 11, 12, 13]

# Gather data
for _, prefix, bits, err, q_curves, _ in techniques:
    for b in n_bits:
        f_path = dataset + '_' + prefix + '_' + str(b) + '.txt'

        with open(f_path) as txt:
            l1 = txt.readline()
            l2 = txt.readline()

            curve = [float(v) for v in l1.split()]
            bbp = np.sum(curve)
            e = float(l2)

            bits.append(bbp)
            err.append(e)
            q_curves.append(curve)

fig, ax = plt.subplots(2, 2)

x = np.arange(len(n_bits))

width = 0.3
offsets = [x - width, x, x + width]

lw = .1

for (tech_name, prefix, bits, err, q_curves, color), offset in zip(techniques, offsets):
    ax[0, 0].bar(offset, err, width, label=tech_name)
    ax[0, 1].bar(offset, bits, width, label=tech_name)

    ax[1, 0].scatter(bits, err, label=tech_name)

    for i, bits in zip(range(len(n_bits)), n_bits):
        curve = q_curves[i]

        if prefix == 'u':
            curve = np.array(q_curves[i]) + lw
        elif prefix == 'utb':
            curve = np.array(q_curves[i])
        elif prefix == 'up':
            curve = np.array(q_curves[i]) - lw

        if prefix != 'up':
            idx = np.arange(len(curve[1:]))
            ax[1, 1].plot(idx, curve[1:], color=color)
        else:
            idx = np.arange(len(curve[1:-1]))
            ax[1, 1].plot(idx, curve[1:-1], color=color)


ax[0, 0].set_yscale('log')
ax[0, 0].set_ylabel('Error')
ax[0, 0].set_xlabel('Starting quantization bits')

ax[0, 1].set_ylabel('Bits per pixel')
ax[0, 1].set_xlabel('Starting quantization bits')

ax[1, 0].set_yscale('log')
ax[1, 0].set_xlabel('Bits per pixel')
ax[1, 0].set_ylabel('Error')

ax[1, 1].set_xlabel('AC component order')
ax[1, 1].set_ylabel('Quantization bits')

ax[0, 0].set_xticks(x, labels=[str(b) for b in n_bits])
ax[0, 1].set_xticks(x, labels=[str(b) for b in n_bits])

ax[0, 0].legend()
ax[0, 1].legend()

fig.tight_layout()


if use_latex:
    plt.savefig('supplemental/plots/chart_and_stuffed_toy_ms.pgf')
else:
    plt.show()
