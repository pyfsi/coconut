from coconut.examples.post_processing.post_processing import *

import numpy as np

# case to be plotted
case_path = 'case_results.pickle'
fsi_case = 'fsi1'

# create PostProcess object
pp = PostProcess(case_path)

# create SubSet of correct point(s)
sx = pp.add_subset(interface='interface_x', model_part='beamrightoutside_nodes')
mask = (abs(sx.get_all_initial_coordinates()[:, 1] - 0.2) < 1e-16)
sx.select_points(mask)

# get displacement values
ux = np.mean(sx.get_values('displacement', component='x'), axis=1)  # y displacement of tip
uy = np.mean(sx.get_values('displacement', component='y'), axis=1)  # x displacement of tip

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

time = sx.get_all_times()

# benchmark data
turek_benchmark = {
    'fsi1':
        {
            'ux':
                {'mean': 2.2700e-5, 'amplitude': 0, 'frequency': 0},
            'uy':
                {'mean': 8.2090e-4, 'amplitude': 0, 'frequency': 0},
            'drag':
                {'mean': 1.4295e+1, 'amplitude': 0, 'frequency': 0},
            'lift':
                {'mean': 7.6380e-1, 'amplitude': 0, 'frequency': 0}
        },
}


coconut = {u_name: {'mean': u[1], 'amplitude': 0, 'frequency': 0} for u_name, u in
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
