import matplotlib.pyplot as plt
import numpy as np
from coconut.examples.post_processing.post_processing import *

# Not excluded in BC
pp = PostProcess('../Alexiades_not_BC_excluded/case_results.pickle')
print(pp)
print("\n")

sx = pp.add_subset(interface='interface_x', model_part='boundary_in_nodes')
sy = pp.add_subset(interface='interface_y', model_part='boundary_out_faces')

# Excluded in BC
pp_excl = PostProcess('../Alexiades_solid_py/case_results.pickle')
print(pp_excl)
print("\n")

sx_excl = pp_excl.add_subset(interface='interface_x', model_part='boundary_in_nodes')
sy_excl = pp_excl.add_subset(interface='interface_y', model_part='boundary_out_faces')

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
    ani_front_excl = Animation2d(sx_excl, 'coordinates', 'coordinates', 'x', 'y', aspect='auto', func_animation_settings=dict(skip=49))
    #ani_front_excl.save('figs_ani/ani_front_excl.gif')

    Animation2d(sy_excl, 'initial_coordinates', 'heat_flux', x_component='y', func_animation_settings=dict(skip=49))
    plt.show()
    plt.close()

# Compare interface at t = 400 s
t_excl = 40000
t_not = 40000
x_not = sx.get_values('coordinates', 'x')
y_not = sx.get_values('coordinates', 'y')
x_not = x_not[t_not,:].flatten()
y_not = y_not[t_not,:].flatten()

x_excl = sx_excl.get_values('coordinates', 'x')
y_excl = sx_excl.get_values('coordinates', 'y')
x_excl = x_excl[t_excl,:].flatten()
y_excl = y_excl[t_excl,:].flatten()

line_not, = plt.plot(x_not, y_not, '-b', label='not excl')
line_excl, = plt.plot(x_excl, y_excl, '-g', label='excluded')
plt.ylabel('y-coordinate [m]')
plt.xlabel('x-coordinate [m]')
#plt.xlim((0, 0.1))
#plt.title('Comparison between BL mesh and non-BL mesh at 10 s')
plt.legend(handles=[line_not, line_excl])
plt.savefig('figs_ani/comp_coupling_not.png')
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