import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from os.path import join, dirname

matplotlib.rcParams['font.size'] = 20
matplotlib.rcParams['axes.linewidth'] = 1
matplotlib.rcParams['lines.linewidth'] = 2
matplotlib.rcParams['axes.labelweight'] = 'bold'
matplotlib.rcParams['axes.labelpad'] = 15
matplotlib.rcParams['xtick.major.size'] = 6
matplotlib.rcParams['xtick.minor.size'] = 2.4
matplotlib.rcParams['ytick.major.size'] = 6
matplotlib.rcParams['ytick.minor.size'] = 2.4
matplotlib.rcParams['savefig.dpi'] = 300
matplotlib.rcParams['figure.figsize'] = 20, 8

w_dir = dirname(__file__)
data = np.loadtxt(join(w_dir, '../CSM/forces_deformed_geometry/SF_deformed_yarn_timestep2000_nodes.txt'), skiprows=1)

data_copy = np.zeros(data.shape)

n_points = 990
for i in range(n_points):
    data_copy[2 * i, :] = data[i, :]
    if i < (n_points - 1):
        data_copy[2 * i + 1, :] = data[n_points + i, :]

x0, y0, z0 = data_copy[:, 1], data_copy[:, 2], data_copy[:, 3]
x1, y1, z1 = data_copy[:, 4] + x0, data_copy[:, 5] + y0, data_copy[:, 6] + z0
print(np.all(x0[:-1] <= x0[1:]))
print(x0, x1)

l0 = [0]
l1 = [0]

for i in range(x0.shape[0] - 1):
    l0.append(l0[i] + np.sqrt((x0[i] - x0[i + 1]) ** 2 + (y0[i] - y0[i + 1]) ** 2 + (z0[i] - z0[i + 1]) ** 2))
    l1.append(l1[i] + np.sqrt((x1[i] - x1[i + 1]) ** 2 + (y1[i] - y1[i + 1]) ** 2 + (z1[i] - z1[i + 1]) ** 2))

x_tube0 = -0.0352
x_air0 = 0.0352
x_air1 = 0.04045
x_tube1 = 0.283102382
x_flow0 = -0.05
x_flow1 = 0.5

i_tube0 = np.argmin(np.abs(x1 - x_tube0))
i_tube1 = np.argmin(np.abs(x1 - x_tube1))
i_air0 = np.argmin(np.abs(x1 - x_air0))
i_air1 = np.argmin(np.abs(x1 - x_air1))
i_flow0 = np.argmin(np.abs(x1 - x_flow0))
i_flow1 = np.argmin(np.abs(x1 - x_flow1))

f_min = np.min(data_copy[:, 7])
f_max = np.max(data_copy[:, 7])

fig, ax_f = plt.subplots()
ax_f.plot(l1, data_copy[:, 7], c='k')

# plot contact body envelope
air_inlet = patches.Rectangle((l0[i_air0], -1), l0[i_air1] - l0[i_air0], 2, alpha=0.5, fc='tab:red', ec='tab:red',
                              fill=True, label='Air inlet')
nozzle = patches.Rectangle((l0[i_tube0], -1), l0[i_tube1] - l0[i_tube0], 2, alpha=0.3, fc='tab:grey', ec='tab:grey',
                           fill=True, label='Contact body')
ax_f.add_patch(nozzle)
ax_f.add_patch(air_inlet)

# plot flow domain boundaries
ax_f.plot([l0[i_flow0], l0[i_flow0]], [-1, 1], c='b', ls=':', lw=2)
ax_f.plot([l0[i_flow1], l0[i_flow1]], [-1, 1], c='b', ls=':', lw=2, label='Flow domain boundary')

ax_f.set_xlabel('Distance along (deformed) yarn centerline [m]')
ax_f.set_ylabel('Axial force [N]', color='k')
plt.tight_layout()
ax_f.set_ylim(2 * f_min, 1.1 * f_max)
ax_f.set_xlim(0.5, 2)
ax_f.legend()
plt.tight_layout()
plt.show()
