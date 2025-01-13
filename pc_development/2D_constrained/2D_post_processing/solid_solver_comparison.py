import matplotlib.pyplot as plt
import numpy as np
from coconut.examples.post_processing.post_processing import *

# With Fluent solid solver
pp_fl = PostProcess('../Alexiades/case_results.pickle')
print(pp_fl)
print("\n")

sx_fl = pp_fl.add_subset(interface='interface_x', model_part='boundary_in_nodes')
sy_fl = pp_fl.add_subset(interface='interface_y', model_part='boundary_out_faces')

# With python_solid_solver
pp_py = PostProcess('../Alexiades_solid_py/case_results.pickle')
print(pp_py)
print("\n")

sx_py = pp_py.add_subset(interface='interface_x', model_part='boundary_in_nodes')
sy_py = pp_py.add_subset(interface='interface_y', model_part='boundary_out_faces')

#print(sx_py.get_all_initial_coordinates())
#print(sx_py.get_values('displacement', 'x'))
#print(sy_py.get_values('heat_flux'))

animations = False
if animations:
    # Animate interface displacement & heat flux of non-conservative simulation
    ani_front = Animation2d(sx_fl, 'coordinates', 'coordinates', 'x', 'y', aspect='auto', func_animation_settings=dict(skip=19))
    #ani_front.save('figs_ani/ani_front_noncons.gif')

    Animation2d(sy_fl, 'initial_coordinates', 'heat_flux', x_component='y', func_animation_settings=dict(skip=19))
    plt.show()
    plt.close()

    # Animate interface displacement & heat flux of conservative simulation
    ani_front_py = Animation2d(sx_py, 'coordinates', 'coordinates', 'x', 'y', aspect='auto', func_animation_settings=dict(skip=19))
    #ani_front_py.save('figs_ani/ani_front_py.gif')

    Animation2d(sy_py, 'initial_coordinates', 'heat_flux', x_component='y', func_animation_settings=dict(skip=19))
    plt.show()
    plt.close()

# Compare interface at t = 100 s
t = 8000
x_fl = sx_fl.get_values('coordinates', 'x')
y_fl = sx_fl.get_values('coordinates', 'y')
x_fl = x_fl[t,:].flatten()
y_fl = y_fl[t,:].flatten()

x_py = sx_py.get_values('coordinates', 'x')
y_py = sx_py.get_values('coordinates', 'y')
x_py = x_py[t,:].flatten()
y_py = y_py[t,:].flatten()

line_fl, = plt.plot(x_fl, y_fl, '-b', label='Solid: Fluent')
line_py, = plt.plot(x_py, y_py, '-g', label='Solid: Python')
plt.ylabel('y-coordinate [m]')
plt.xlabel('x-coordinate [m]')
plt.xlim((0, 0.1))
#plt.title('Comparison between BL mesh and non-BL mesh at 10 s')
plt.legend(handles=[line_fl, line_py])
plt.savefig('figs_ani/comp_py_fluent_solid_solver.png')
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