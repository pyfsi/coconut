from coconut.examples.post_processing.post_processing import *

import numpy as np
import matplotlib.pyplot as plt

# case to be plotted
case_path = 'case_results.pickle'
fsi_case = 'fsi2'

# create PostProcess object
pp = PostProcess(case_path)

# create SubSet of correct point(s)
sx = pp.add_subset(interface='interface_x', model_part='beamrightoutside_nodes')
mask = (abs(sx.get_all_initial_coordinates()[:, 1] - 0.2) < 1e-16)
sx.select_points(mask)

# get displacement values
ux = np.mean(sx.get_values('displacement', component='x'), axis=1)  # y displacement of tip
uy = np.mean(sx.get_values('displacement', component='y'), axis=1)  # x displacement of tip

time = sx.get_all_times()

drag = []
lift = []
with open(f'CFD/forces.frp', 'r') as f:
    for line in f:
        if line.startswith('Net'):
            forces = line.split()
            if forces[1].startswith('('):
                drag.append(float(forces[7].strip('()')))
                lift.append(float(forces[8].strip('()')))
if len(drag) > time.size:
    drag = drag[:time.size]
    lift = lift[:time.size]
elif len(drag) < time.size:
    drag = [np.nan] * (time.size - len(drag)) + drag
    lift = [np.nan] * (time.size - len(lift)) + lift
drag = np.array(drag) * 100  # mesh depth (z-direction) is 0.01 m
lift = np.array(lift) * 100  # mesh depth (z-direction) is 0.01 m

# benchmark data
turek_benchmark = {
    'fsi2':
        {
            'ux':
                {'mean': -1.4580e-2, 'amplitude': 1.2440e-2, 'frequency': 3.8},
            'uy':
                {'mean': 1.2300e-3, 'amplitude': 8.0600e-2, 'frequency': 2.0},
            'drag':
                {'mean': 2.0883e+2, 'amplitude': 7.3750e+1, 'frequency': 3.8},
            'lift':
                {'mean': 8.8000e-1, 'amplitude': 2.3420e+2, 'frequency': 2.0},
        }
}


# calculate coconut data
def determine_data(u, time):
    if np.all(u == np.zeros_like(time)):
        return {'mean': 0, 'amplitude': 0, 'frequency': 0}
    else:
        n = min(int(0.12 * time.size), sum(~np.isnan(u)))
        u_lim = u[-n:]  # only look at last part
        time_lim = time[-n:]
        u_m = np.mean(u_lim)
        u_switch = np.where(np.diff(u_lim < u_m) != 0)[0]  # determine indices of mean crossings
        if u_switch.size < 3:
            raise RuntimeError('Not enough time step data to determine the mean, amplitude and frequency')
        u_per = u_lim[u_switch[-3]:u_switch[-1] + 1]  # last period of u (between two mean crossings)
        time_per = time_lim[u_switch[-3]:u_switch[-1] + 1]
        u_max = np.max(u_per)
        u_min = np.min(u_per)
        u_mean = 0.5 * (u_max + u_min)
        u_ampl = 0.5 * (u_max - u_min)
        u_frequency = 1 / (time_per[-1] - time_per[0])
    return {'mean': u_mean, 'amplitude': u_ampl, 'frequency': u_frequency}


stop = 35  # include simulation time up to 'stop' seconds
coconut = {u_name: determine_data(u[time <= stop], time[time <= stop]) for u_name, u in
           (('ux', ux), ('uy', uy), ('drag', drag), ('lift', lift))}

# summary
for var1_str, var1_name, var2_str, var2_name in (('ux', 'displacement x', 'uy', 'displacement y'),
                                                 ('drag', 'drag', 'lift', 'lift')):
    out = f'{fsi_case:12s}{var1_name:34s}{var2_name:34s}\n' + '-' * 80 + '\n'
    for dict_name, dct in (('benchmark', turek_benchmark[fsi_case]), ('coconut', coconut)):
        out += f'{dict_name:12s}'
        for var in (var1_str, var2_str):
            out += f"{dct[var]['mean']:11.4e} +/-{dct[var]['amplitude']:.4e} [{dct[var]['frequency']:.1f}]"
            out += '  ' if var == var1_str else '\n'
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
