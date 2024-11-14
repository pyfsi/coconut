import matplotlib.pyplot as plt
import numpy as np
from coconut.examples.post_processing.post_processing import *

# With non-conservative mapper
pp = PostProcess('case_results/case_non_conservative.pickle')
print(pp)
print("\n")

sx = pp.add_subset(interface='interface_x', model_part='boundary_in_nodes')
sy = pp.add_subset(interface='interface_y', model_part='boundary_out_faces')

# With conservative mapper
pp_c = PostProcess('case_results/case_conservative.pickle')
print(pp_c)
print("\n")

sx_c = pp_c.add_subset(interface='interface_x', model_part='boundary_in_nodes')
sy_c = pp_c.add_subset(interface='interface_y', model_part='boundary_out_faces')

#print(sx.get_all_initial_coordinates())
#print(sx.get_values('displacement', 'x'))
#print(sy.get_values('heat_flux'))

# Animate interface displacement & heat flux of non-conservative simulation
ani_front = Animation2d(sx, 'coordinates', 'coordinates', 'x', 'y', aspect='auto', func_animation_settings=dict(skip=299))
#ani_front.save('figs_ani/ani_front_noncons.gif')

Animation2d(sy, 'initial_coordinates', 'heat_flux', x_component='y', func_animation_settings=dict(skip=299))
plt.show()
plt.close()

# Animate interface displacement & heat flux of conservative simulation
ani_front_cons = Animation2d(sx_c, 'coordinates', 'coordinates', 'x', 'y', aspect='auto', func_animation_settings=dict(skip=299))
#ani_front_cons.save('figs_ani/ani_front_cons.gif')

Animation2d(sy_c, 'initial_coordinates', 'heat_flux', x_component='y', func_animation_settings=dict(skip=299))
plt.show()
plt.close()

# Compare interface at t = 480 s
x_nc = sx.get_values('coordinates', 'x')
y_nc = sx.get_values('coordinates', 'y')
x_nc = x_nc[48001,:].flatten()
y_nc = y_nc[48001,:].flatten()

x_c = sx_c.get_values('coordinates', 'x')
y_c = sx_c.get_values('coordinates', 'y')
x_c = x_c[48001,:].flatten()
y_c = y_c[48001,:].flatten()

line_nc, = plt.plot(x_nc, y_nc, '-b', label='non-conservative')
line_c, = plt.plot(x_c, y_c, '-g', label='conservative')
plt.ylabel('y-coordinate [m]')
plt.xlabel('x-coordinate [m]')
#plt.xlim((0, 0.1))
#plt.title('Comparison between n2f mappers at 480 s')
plt.legend(handles=[line_nc, line_c])
plt.savefig('figs_ani/cons_comp_zoomed.png')
plt.show()
plt.close()


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