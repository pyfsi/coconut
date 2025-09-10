import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy import integrate

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

data_coords = np.loadtxt('../CFD/coordinates_update_timestep2000.dat', skiprows=1)

l1 = [0]

for i in range(data_coords.shape[0] - 1):
    l1.append(l1[i] + np.sqrt(
        (data_coords[i, 0] - data_coords[i + 1, 0]) ** 2 + (data_coords[i, 1] - data_coords[i + 1, 1]) ** 2 + (
                    data_coords[i, 2] - data_coords[i + 1, 2]) ** 2))

data_forces = np.loadtxt('../CFD/traction_timestep2000.dat', skiprows=1)[:, :-1]

tangents = data_coords[:, 3:]
tangents /= np.linalg.norm(tangents, axis=1).reshape(-1, 1)
print(data_forces.shape, tangents.shape)

axial_force = np.zeros(tangents.shape[0])
for i in range(tangents.shape[0]):
    axial_force[i] = np.dot(data_forces[i, :], tangents[i, :])
total_force = np.linalg.norm(data_forces, axis=1)

print(axial_force.shape, len(l1))

x_tube0 = -0.0352
x_air0 = 0.0352
x_air1 = 0.04045
x_tube1 = 0.283102382
x_flow0 = -0.05
x_flow1 = 0.5

i_tube0 = np.argmin(np.abs(data_coords[:, 0] - x_tube0))
i_tube1 = np.argmin(np.abs(data_coords[:, 0] - x_tube1))
i_air0 = np.argmin(np.abs(data_coords[:, 0] - x_air0))
i_air1 = np.argmin(np.abs(data_coords[:, 0] - x_air1))
i_flow0 = np.argmin(np.abs(data_coords[:, 0] - x_flow0))
i_flow1 = np.argmin(np.abs(data_coords[:, 0] - x_flow1))

force_ALM = integrate.trapezoid(data_forces[i_flow0:i_flow1, :], np.array(l1[i_flow0:i_flow1]).reshape(-1, 1), axis=0)
print('Line integral of forces: ', force_ALM)

f_min = np.min(axial_force)
f_max = np.max(axial_force)
print(f_min, f_max)

plt.figure()
plt.plot(l1, total_force, c='k', label='Total force (magnitude)')
plt.plot(l1, axial_force, c='g', label='Axial force')

# plot contact body envelope
ax = plt.gca()
air_inlet = patches.Rectangle((l1[i_air0], f_min-1), l1[i_air1] - l1[i_air0], f_max+abs(f_min)+2, alpha=0.5, fc='tab:red', ec='tab:red',
                              fill=True, label='Air inlet')
nozzle = patches.Rectangle((l1[i_tube0], f_min-1), l1[i_tube1] - l1[i_tube0], f_max+abs(f_min)+2, alpha=0.3, fc='tab:grey', ec='tab:grey',
                           fill=True, label='Contact body')
ax.add_patch(nozzle)
ax.add_patch(air_inlet)

# plot flow domain boundaries
plt.plot([l1[i_flow0], l1[i_flow0]], [f_min-1, f_max+1], c='b', ls=':', lw=2)
plt.plot([l1[i_flow1], l1[i_flow1]], [f_min-1, f_max+1], c='b', ls=':', lw=2, label='Flow domain boundary')

plt.xlabel('Distance along (deformed) yarn centerline [m]')
plt.ylabel('Aerodynamic force [N/m]')
plt.tight_layout()
plt.ylim(2 * f_min, 1.1 * f_max)
plt.xlim(0.5, 2)
plt.legend()
plt.tight_layout()
plt.show()
