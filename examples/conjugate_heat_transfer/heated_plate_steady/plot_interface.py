import matplotlib.pyplot as plt
from coconut.examples.post_processing.post_processing import *

pp = PostProcess('case_results.pickle')
print(pp)

sx = pp.add_subset(model_part='boundary_in_faces')
sy = pp.add_subset(model_part='boundary_out_faces')

Plot2d(sx, 'initial_coordinates', 'temperature', x_component='x', time_step=1)
plt.show()

Plot2d(sy, 'initial_coordinates', 'heat_flux', x_component='x', time_step=1)
plt.show()