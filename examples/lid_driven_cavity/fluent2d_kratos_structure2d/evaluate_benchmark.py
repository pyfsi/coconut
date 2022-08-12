from coconut.examples.post_processing.animate import AnimationFigure

import numpy as np
import matplotlib.pyplot as plt
import pickle
import os

# different cases to be plotted
common_path = ''
case_paths = ['case_results.pickle']
legend_entries = ['CoCoNuT']
data = {}

for name, path in zip(legend_entries, case_paths):
    # case to be plotted
    case_path = os.path.join(common_path, path)

    # load case
    results = pickle.load(open(case_path, 'rb'))

    # make figure and create animation for each case
    animation_figure = AnimationFigure()  # figure for coordinate animations

    solution = results['solution_x']
    interface = results['interface_x']
    dt = results['delta_t']
    time_step_start = results['timestep_start']
    # create animate object
    try:
        animation = animation_figure.add_animation(solution, interface, dt, time_step_start, variable='displacement',
                                                   model_part_name='cavitybottom_nodes')
    except Exception:
        animation = animation_figure.add_animation(solution, interface, dt, time_step_start, variable='displacement',
                                                   model_part_name='StructureInterface2D_StructureInterface_output')

    # select points and component of variable to plot
    coordinates = animation.coordinates
    mask_x = (abs(coordinates[:, 0] - 0.5) < 1e-16)
    mask_y = (coordinates[:, 1] > -np.inf)
    mask_z = (coordinates[:, 2] > -np.inf)
    abscissa = 0  # x-axis
    component = 1  # y-component

    animation.initialize(mask_x, mask_y, mask_z, abscissa, component)

    # variables
    uy = np.mean(np.array(animation.solution), axis=1)  # y displacement of center
    ux = np.mean(np.array(animation.displacement), axis=1)  # x displacement of center
    time = animation.time_step_start + np.arange(animation.time_steps + 1) * animation.dt

    data.update({name: (time, ux, uy)})

# plot
save = False
xlim = (0, 70)
ylim = (-0.05, 0.30)
colors = ['blue', 'red', 'green', 'orange', 'purple']

_, ax = plt.subplots(figsize=(10, 7))
j = 0
time = np.linspace(0, 70, 701)
for reference in ('Mok', 'Valdes'):
    time, uy = np.loadtxt(f'../{reference.lower()}.csv', skiprows=1, delimiter=',', unpack=True)
    plt.plot(time, uy, label=reference, color=colors[j], linewidth=1.5, linestyle='-', marker='o', markersize=3)
    j += 1

for i, name in enumerate(legend_entries):
    time, ux, uy = data[name]
    # plt.plot(time, ux, label=name + ' displacement x', color=colors[j + i], linewidth=1.5, linestyle='--')
    plt.plot(time, uy, label=name + ' displacement y', color=colors[j + i], linewidth=1.5, linestyle='-')

plt.xlabel('time (s)')
plt.ylabel('y-displacement (m)')
plt.xlim(*xlim)
plt.ylim(*ylim)
ax.tick_params(axis='both', direction='in', pad=8, top=True, right=True)
plt.tight_layout()
plt.legend(loc='lower right')
if save:
    plt.savefig('comparison_fluent.png', dpi=300)

plt.show()
