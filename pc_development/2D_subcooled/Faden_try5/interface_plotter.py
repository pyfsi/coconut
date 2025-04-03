import matplotlib.pyplot as plt
import numpy as np
from coconut.examples.post_processing.post_processing import *

# Reference
pp = PostProcess('case_results.pickle')
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
    ani_front_new = Animation2d(sx_new, 'coordinates', 'coordinates', 'x', 'y', aspect='auto', func_animation_settings=dict(skip=19))
    #ani_front_new.save('figs_ani/ani_front_new_ts.gif')

    Animation2d(sy_new, 'initial_coordinates', 'heat_flux', x_component='y', func_animation_settings=dict(skip=19))
    plt.show()
    plt.close()

    # Animate interface displacement & heat flux of conservative simulation
    ani_front = Animation2d(sx, 'coordinates', 'coordinates', 'x', 'y', aspect='auto', func_animation_settings=dict(skip=19))
    #ani_front.save('figs_ani/ani_front.gif')

    Animation2d(sy, 'initial_coordinates', 'heat_flux', x_component='y', func_animation_settings=dict(skip=19))
    plt.show()
    plt.close()

# Compare interface
t = 28.8 # s
dt = 0.1 # s

x = sx.get_values('coordinates', 'x')
y = sx.get_values('coordinates', 'y')
x = x[int(t/dt),:].flatten()
y = y[int(t/dt),:].flatten()

line, = plt.plot(x, y, '-k', label='Reference')
plt.ylabel('y-coordinate [m]')
plt.xlabel('x-coordinate [m]')
plt.xlim((0, 0.04))
plt.legend(handles=[line])
#plt.savefig('interface.png')
plt.show()
plt.close()