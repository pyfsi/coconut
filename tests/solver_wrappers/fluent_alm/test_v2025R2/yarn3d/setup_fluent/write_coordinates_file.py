import numpy as np

x = np.linspace(start=-0.25 - 1e-6, stop=0.25 + 1e-6, num=12)
coordinates = np.zeros((x.shape[0], 3))
coordinates[:, 0] = x

np.savetxt('coordinates_timestep0.dat', coordinates, fmt='%27.17e%27.17e%27.17e', header=f'{coordinates.shape[0]}',
           comments='')
