from coconut.examples.post_processing.animate import AnimationFigure

import numpy as np
import pickle

# case to be plotted
case_path = 'case_results.pickle'
fsi_case = 'fsi1'

# load case
results = pickle.load(open(case_path, 'rb'))

# make figure and create animation for each case
animation_figure = AnimationFigure()  # figure for coordinate animations

solution = results['solution_x']
interface = results['interface_x']
dt = results['delta_t']
time_step_start = results['timestep_start']
# create animate object
animation = animation_figure.add_animation(solution, interface, dt, time_step_start, variable='displacement',
                                           name=fsi_case, model_part_name='beamrightoutside_nodes')
# select points and component of variable to plot
coordinates = animation.coordinates
mask_x = (coordinates[:, 0] > -np.inf)
mask_y = (abs(coordinates[:, 1] - 0.2) < 1e-16)
mask_z = (coordinates[:, 2] > -np.inf)
abscissa = 0  # x-axis, not used
component = 1  # y-component, not used

animation.initialize(mask_x, mask_y, mask_z, abscissa, component)

# variables
uy = np.mean(np.array(animation.solution), axis=1)  # y displacement of tip
ux = np.mean(np.array(animation.displacement), axis=1)  # x displacement of tip
drag = []
lift = []
with open(f'CFD/forces.frp', 'r') as f:
    for line in f:
        if line.startswith('Net'):
            forces = line.split()
            if forces[1].startswith('('):
                drag.append(float(forces[7].strip('()')))
                lift.append(float(forces[8].strip('()')))
drag = np.array(drag)
lift = np.array(lift)
time = animation.time_step_start + np.arange(animation.time_steps + 1) * animation.dt

# benchmark data
turek_benchmark = {
    'fsi1':
        {
            'ux':
                {'mean': 2.2700e-5, 'amplitude': 0, 'period': 0},
            'uy':
                {'mean': 8.2090e-4, 'amplitude': 0, 'period': 0},
            'drag':
                {'mean': 1.4295e+1, 'amplitude': 0, 'period': 0},
            'lift':
                {'mean': 7.6380e-1, 'amplitude': 0, 'period': 0}
        },
}


coconut = {u_name: {'mean': u[1], 'amplitude': 0, 'period': 0} for u_name, u in
           (('ux', ux), ('uy', uy), ('drag', drag), ('lift', lift))}

# summary
for var1_str, var1_name, var2_str, var2_name in (('ux', 'displacement x', 'uy', 'displacement y'),
                                                 ('drag', 'drag', 'lift', 'lift')):
    out = f"{fsi_case}\t\t{var1_name:30s}\t\t{var2_name:30s}\n" + '-' * 80 + '\n'
    for dict_name, dict in (('benchmark', turek_benchmark[fsi_case]), ('coconut', coconut)):
        out += f"{dict_name:12}"
        for var in (var1_str, var2_str):
            out += f"{dict[var]['mean']:11.4e} +/-{dict[var]['amplitude']:.4e} [{dict[var]['period']:.1f}]"
            out += '\t\t' if var == var1_str else '\n'
    print(out)
