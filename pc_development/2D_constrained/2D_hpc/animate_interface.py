import matplotlib.pyplot as plt
import numpy as np
from coconut.examples.post_processing.post_processing import *

pp = PostProcess('case_results.pickle')
print(pp)
print("\n")
#pp.print_summary()

sx = pp.add_subset(interface='interface_x', model_part='boundary_in_nodes')
sy = pp.add_subset(interface='interface_y', model_part='boundary_out_faces')

#print(sx.get_all_initial_coordinates())
#print(sx.get_values('displacement', 'x'))
#print(sy.get_values('heat_flux'))

# Animate interface displacement & heat flux
ani_front = Animation2d(sx, 'coordinates', 'coordinates', 'x', 'y', aspect='auto', func_animation_settings=dict(skip=9))
ani_plot = Plot2d(sx, 'coordinates', 'coordinates', 'x', 'y', aspect='auto', time_step=580)
#ani_front.save('figs_ani/ani_front.gif')
plt.show()
plt.close()
#ani_disp = Animation2dDisplacement(sx, x_component='y', y_component='x')

Animation2d(sy, 'initial_coordinates', 'heat_flux', x_component='y', func_animation_settings=dict(skip=9))
plt.show()

"""
# Time evolution of average temperature & heat flux
temp = sx.get_values('temperature')[:,:,0]
hf = sy.get_values('heat_flux')[:,:,0]
temp_avg = np.mean(temp, axis=1)
hf_avg = np.mean(hf, axis=1)

dt = pp.get_data()['delta_t']
ts_start = pp.get_data()['timestep_start']
nr_ts = np.shape(temp)[0]
time = np.array(range(ts_start, nr_ts, 1))*dt

plt.plot(time, temp_avg)
plt.ylabel('Average temperature [K]')
plt.xlabel('Time [s]')
plt.show()

plt.plot(time, -hf_avg)
plt.ylabel('Average heat flux [W/m^2]')
plt.xlabel('Time [s]')
plt.show()
"""