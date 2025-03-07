import matplotlib.pyplot as plt
import numpy as np
from coconut.examples.post_processing.post_processing import *

# With explicit coupling
pp = PostProcess('../Alexiades_explicit/case_results.pickle')
print(pp)
print("\n")

sx = pp.add_subset(interface='interface_x', model_part='boundary_in_nodes')
sy = pp.add_subset(interface='interface_y', model_part='boundary_out_faces')

# With aitken coupling
pp_c = PostProcess('../Alexiades_solid_py/case_results.pickle')
print(pp_c)
print("\n")

sx_c = pp_c.add_subset(interface='interface_x', model_part='boundary_in_nodes')
sy_c = pp_c.add_subset(interface='interface_y', model_part='boundary_out_faces')

#print(sx.get_all_initial_coordinates())
#print(sx.get_values('displacement', 'x'))
#print(sy.get_values('heat_flux'))

animations = False
if animations:
    # Animate interface displacement & heat flux of non-conservative simulation
    ani_front = Animation2d(sx, 'coordinates', 'coordinates', 'x', 'y', aspect='auto', func_animation_settings=dict(skip=49))
    #ani_front.save('figs_ani/ani_front_noncons.gif')

    Animation2d(sy, 'initial_coordinates', 'heat_flux', x_component='y', func_animation_settings=dict(skip=49))
    plt.show()
    plt.close()

    # Animate interface displacement & heat flux of conservative simulation
    ani_front_c = Animation2d(sx_c, 'coordinates', 'coordinates', 'x', 'y', aspect='auto', func_animation_settings=dict(skip=49))
    #ani_front_c.save('figs_ani/ani_front_c.gif')

    Animation2d(sy_c, 'initial_coordinates', 'heat_flux', x_component='y', func_animation_settings=dict(skip=49))
    plt.show()
    plt.close()

# Compare interface at given time
time = [100, 200, 450, 700]
lines = []

x_expl = sx.get_values('coordinates', 'x')
y_expl = sx.get_values('coordinates', 'y')
x_c = sx_c.get_values('coordinates', 'x')
y_c = sx_c.get_values('coordinates', 'y')

for t in time:
    x_expl_tmp = x_expl[t*100,:].flatten()
    y_expl_tmp = y_expl[t*100,:].flatten()
    x_c_tmp = x_c[t*100,:].flatten()
    y_c_tmp = y_c[t*100,:].flatten()

    line_expl, = plt.plot(x_expl_tmp, y_expl_tmp, '--r', label='explicit')
    lines.append(line_expl)
    line_c, = plt.plot(x_c_tmp, y_c_tmp, '-k', label='aitken')
    lines.append(line_c)

legend_handles = [plt.plot([], [], '--r')[0], plt.plot([], [], '-k')[0]] # Get handles
legend_labels = ['explicit', 'Aitken']

plt.ylabel('y-coordinate [m]')
plt.xlabel('x-coordinate [m]')
plt.xlim((0, 0.05))
plt.legend(legend_handles, legend_labels)
plt.tight_layout()
plt.savefig('figs_ani/comp_coupling_expl.png')
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