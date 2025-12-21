import matplotlib.pyplot as plt
import numpy as np
import io
from coconut.examples.post_processing.post_processing import *

file = 'thick_overset_ramped/case_results.pickle'
time_step = 100

pp = PostProcess(file)
sx = pp.add_subset(interface='interface_x', model_part='inner_in_nodes')
sy = pp.add_subset(interface='interface_y', model_part='inner_out_faces')

x = sx.get_values('coordinates', 'x')
y = sx.get_values('coordinates', 'y')
x = x[time_step, :].flatten()
y = y[time_step, :].flatten()

plt.figure()
plt.plot(x * 1000, y * 1000, 'b.')
plt.xlabel('x-coordinate (mm)')
plt.ylabel('y-coordinate (mm)')
plt.grid(True)
plt.show()