import matplotlib.pyplot as plt
import numpy as np
from coconut.examples.post_processing.post_processing import *

# With large time steps
pp_large = PostProcess('../Alexiades_large_ts/case_results.pickle')
print(pp_large)
print("\n")

sx_large = pp_large.add_subset(interface='interface_x', model_part='boundary_in_nodes')
sy_large = pp_large.add_subset(interface='interface_y', model_part='boundary_out_faces')

# With python_solid_solver
pp = PostProcess('../Alexiades_solid_py/case_results.pickle')
print(pp)
print("\n")

sx = pp.add_subset(interface='interface_x', model_part='boundary_in_nodes')
sy = pp.add_subset(interface='interface_y', model_part='boundary_out_faces')

#print(sx.get_all_initial_coordinates())
#print(sx.get_values('displacement', 'x'))
#print(sy.get_values('heat_flux'))

animations = False
if animations:
    # Animate interface displacement & heat flux of non-conservative simulation
    ani_front_large = Animation2d(sx_large, 'coordinates', 'coordinates', 'x', 'y', aspect='auto', func_animation_settings=dict(skip=19))
    #ani_front_large.save('figs_ani/ani_front_large_ts.gif')

    Animation2d(sy_large, 'initial_coordinates', 'heat_flux', x_component='y', func_animation_settings=dict(skip=19))
    plt.show()
    plt.close()

    # Animate interface displacement & heat flux of conservative simulation
    ani_front = Animation2d(sx, 'coordinates', 'coordinates', 'x', 'y', aspect='auto', func_animation_settings=dict(skip=19))
    #ani_front.save('figs_ani/ani_front.gif')

    Animation2d(sy, 'initial_coordinates', 'heat_flux', x_component='y', func_animation_settings=dict(skip=19))
    plt.show()
    plt.close()

# Compare interface at t = 35 s
t = 350
x_large = sx_large.get_values('coordinates', 'x')
y_large = sx_large.get_values('coordinates', 'y')
x_large = x_large[t*10,:].flatten()
y_large = y_large[t*10,:].flatten()

x = sx.get_values('coordinates', 'x')
y = sx.get_values('coordinates', 'y')
x = x[t*100,:].flatten()
y = y[t*100,:].flatten()

line_large, = plt.plot(x_large, y_large, '-b', label='Large TS: 0.1 s')
line, = plt.plot(x, y, '-g', label='Small TS: 0.01 s')
plt.ylabel('y-coordinate [m]')
plt.xlabel('x-coordinate [m]')
plt.xlim((0, 0.1))
#plt.title('Comparison between BL mesh and non-BL mesh at 10 s')
plt.legend(handles=[line_large, line])
plt.savefig('figs_ani/comp_time_step_size.png')
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