from coconut.examples.post_processing.post_processing import *

import numpy as np
import matplotlib.pyplot as plt
import pickle
import os

# different cases to be plotted
common_path = ''
case_paths = ['case_results.pickle']
legend_entries = ['CoCoNuT']
data = {}

pps = [None] * len(case_paths)
sxs = [None] * len(case_paths)
for i, path in enumerate(case_paths):
    pps[i] = pp = PostProcess(path)
    sxs[i] = sx = pp.add_subset(interface='interface_x')

    # select point to plot
    coordinates = sx.get_all_initial_coordinates()
    sx.select_points(abs(coordinates[:, 0] - 0.5) < 1e-16)

    # variables
    ux = np.mean(sx.get_values('displacement', 'x'), axis=1)  # x displacement of center
    uy = np.mean(sx.get_values('displacement', 'y'), axis=1)  # y displacement of center
    time = sx.get_all_times()

    data.update({legend_entries[i]: (time, ux, uy)})

# plot
save = False
xlim = (0, 70)
ylim = (-0.05, 0.30)
colors = ['blue', 'red', 'green', 'orange', 'purple']

_, ax = plt.subplots(figsize=(10, 7))
j = 0
time = np.linspace(0, 70, 701)
for reference in ('Mok', 'Valdes'):
    time, uy = np.loadtxt(f'../{reference.lower()}.csv', skiprows=1, delimiter=',', unpack=True)
    plt.plot(time, uy, label=reference, color=colors[j], linewidth=1.5, linestyle='-', marker='o', markersize=3)
    j += 1

for i, name in enumerate(legend_entries):
    time, ux, uy = data[name]
    plt.plot(time, ux, label=name + ' displacement x', color=colors[j + i], linewidth=1.5, linestyle='--')
    plt.plot(time, uy, label=name + ' displacement y', color=colors[j + i], linewidth=1.5, linestyle='-')

plt.xlabel('time (s)')
plt.ylabel('y-displacement (m)')
plt.xlim(*xlim)
plt.ylim(*ylim)
ax.tick_params(axis='both', direction='in', pad=8, top=True, right=True)
plt.tight_layout()
plt.legend(loc='lower right')
if save:
    plt.savefig('comparison_openfoam.png', dpi=300)

plt.show()
