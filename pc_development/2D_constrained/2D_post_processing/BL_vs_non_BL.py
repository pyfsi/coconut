import matplotlib.pyplot as plt
import numpy as np
from coconut.examples.post_processing.post_processing import *

# With no boundary layer
pp = PostProcess('../Alexiades_only_triangles/case_results.pickle')
print(pp)
print("\n")

sx = pp.add_subset(interface='interface_x', model_part='boundary_in_nodes')
sy = pp.add_subset(interface='interface_y', model_part='boundary_out_faces')

# With boundary layer
pp_BL = PostProcess('../Alexiades_solid_py/case_results.pickle')
print(pp_BL)
print("\n")

sx_BL = pp_BL.add_subset(interface='interface_x', model_part='boundary_in_nodes')
sy_BL = pp_BL.add_subset(interface='interface_y', model_part='boundary_out_faces')

#print(sx.get_all_initial_coordinates())
#print(sx.get_values('displacement', 'x'))
#print(sy.get_values('heat_flux'))

animations = False
if animations:
    # Animate interface displacement & heat flux of non-conservative simulation
    ani_front = Animation2d(sx, 'coordinates', 'coordinates', 'x', 'y', aspect='auto', func_animation_settings=dict(skip=19))
    #ani_front.save('figs_ani/ani_front_noncons.gif')

    Animation2d(sy, 'initial_coordinates', 'heat_flux', x_component='y', func_animation_settings=dict(skip=19))
    plt.show()
    plt.close()

    # Animate interface displacement & heat flux of conservative simulation
    ani_front_BL = Animation2d(sx_BL, 'coordinates', 'coordinates', 'x', 'y', aspect='auto', func_animation_settings=dict(skip=19))
    #ani_front_BL.save('figs_ani/ani_front_BL.gif')

    Animation2d(sy_BL, 'initial_coordinates', 'heat_flux', x_component='y', func_animation_settings=dict(skip=19))
    plt.show()
    plt.close()

# Compare interface at t = 300 s
x_non_BL = sx.get_values('coordinates', 'x')
y_non_BL = sx.get_values('coordinates', 'y')
x_non_BL = x_non_BL[30000,:].flatten()
y_non_BL = y_non_BL[30000,:].flatten()

x_BL = sx_BL.get_values('coordinates', 'x')
y_BL = sx_BL.get_values('coordinates', 'y')
x_BL = x_BL[30000,:].flatten()
y_BL = y_BL[30000,:].flatten()

line_non_BL, = plt.plot(x_non_BL, y_non_BL, '-b', label='no BL')
line_BL, = plt.plot(x_BL, y_BL, '-g', label='BL')
plt.ylabel('y-coordinate [m]')
plt.xlabel('x-coordinate [m]')
plt.xlim((0, 0.1))
#plt.title('Comparison between BL mesh and non-BL mesh at 10 s')
plt.legend(handles=[line_non_BL, line_BL])
plt.savefig('figs_ani/comp_BL_zoom_out.png')
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