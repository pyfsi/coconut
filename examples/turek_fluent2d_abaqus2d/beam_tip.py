from coconut.examples.post_processing.animate import AnimationFigure

import numpy as np
import matplotlib.pyplot as plt
import pickle

# case to be plotted
case_path = 'case_results.pickle'
fsi_case = 'fsi2'

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
drag = [np.nan] * (17501 - len(drag)) + drag
lift = [np.nan] * (17501 - len(lift)) + lift
drag = np.array(drag)
lift = np.array(lift)
time = animation.time_step_start + np.arange(animation.time_steps + 1) * animation.dt

# benchmark data
turek_benchmark = {
    'fsi2':
        {
            'ux':
                {'mean': -1.4580e-2, 'amplitude': 1.2440e-2, 'period': 3.8},
            'uy':
                {'mean': 1.2300e-3, 'amplitude': 8.0600e-2, 'period': 2.0},
            'drag':
                {'mean': 2.0883e+2, 'amplitude': 7.3750e+1, 'period': 3.8},
            'lift':
                {'mean': 8.8000e-1, 'amplitude': 2.3420e+2, 'period': 2.0},
        }
}


# calculate coconut data
def determine_data(u, time):
    if np.all(u == np.zeros_like(time)):
        return {'mean': 0, 'amplitude': 0, 'period': 0}
    else:
        n = min(int(0.12 * time.size), sum(~np.isnan(u)))
        u_lim = u[-n:]  # only look at last part
        time_lim = time[-n:]
        u_m = np.mean(u_lim)
        u_switch = np.where(np.diff(u_lim < u_m) != 0)[0]  # determine indices of mean crossings
        if len(u_switch < 3):
            raise RuntimeError('Not enough time step data to determine the mean, amplitude and period')
        u_per = u_lim[u_switch[-3]:u_switch[-1]]  # last period of u (between two mean crossings)
        time_per = time_lim[u_switch[-3]:u_switch[-1]]
        u_max = np.max(u_per)
        time_u_max = time_per[np.where(u_per == u_max)][0]
        u_min = np.min(u_per)
        time_u_min = time_per[np.where(u_per == u_min)][0]
        u_mean = 0.5 * (u_max + u_min)
        u_ampl = 0.5 * (u_max - u_min)
        u_period = 1 / (2 * abs(time_u_max - time_u_min))
    return {'mean': u_mean, 'amplitude': u_ampl, 'period': u_period}


stop = 34  # include simulation time up to 'stop' seconds
coconut = {u_name: determine_data(u[time <= stop], time[time <= stop]) for u_name, u in
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

# plots
xlim = (stop - 1, stop)  # (34, 35)
save = False
ylim = {'ux': (-0.03, 0), 'uy': (-0.08, 0.1), 'drag': (120, 300), 'lift': (-250, 250)}

for var, var_str, var_name, file_name in ((ux, 'ux', 'displacement x', 'displacement_x.png'),
                                          (uy, 'uy', 'displacement y', 'displacement_y.png'),
                                          (drag, 'drag', 'drag', 'drag.png'),
                                          (lift, 'lift', 'lift', 'lift.png')):
    _, ax = plt.subplots(figsize=(5, 3.5))
    plt.plot(time, var, label='coconut', color='black', linewidth=0.5)
    plt.title(fsi_case.upper())
    plt.xlabel('time')
    plt.ylabel(var_name)
    plt.xlim(*xlim)
    plt.ylim(*ylim[var_str])
    ax.tick_params(axis='both', direction='in', pad=8, top=True, right=True)
    plt.tight_layout()
    plt.legend(loc='upper right')
    if save:
        plt.savefig(file_name, dpi=300)

plt.show()
