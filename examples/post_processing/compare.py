from coconut import tools

import numpy as np
import os
import pickle

# This script compares the average norm of the difference per time step between several cases,
# as wel as the number of coupling iterations in each time step, under the assumption that the same solvers have been
# used. Hence, it is mainly meant to compare the performance of coupling algorithms.
# To generate a result file, include a non-zero int 'write_results' in the settings of the coupled solver.
# Give a name to the case by including the string {'case_name': 'a_name'} in the settings of the coupled solver.

# different cases to be plotted
common_path = '../../examples/tube/'
case_paths = ['tube_flow_tube_structure/case_results.pickle']
legend_entries = ['case_results']

# load cases
results = {}
for name, path in zip(legend_entries, case_paths):
    with open(os.path.join(common_path, path), 'rb') as file:
        results.update({name: pickle.load(file)})

# reference case
case_reference = legend_entries[0]

# check equal number of time steps
time_steps = len(results[case_reference]['iterations'])
for case in legend_entries:
    if not len(results[case]['iterations']) == time_steps:
        raise Exception(f"Number of time steps for case {case} ({len(results[case]['iterations'])}) "
                        f"differs from number of time steps of reference case ({time_steps}).")

for case in legend_entries:
    norm = []
    for i in range(time_steps):
        norm = norm + [np.linalg.norm(results[case_reference]['solution_x'][:, i]
                                      - results[case]['solution_x'][:, i])]
    norm_average = np.array(norm).mean()
    tools.print_info(f'Average norm of the difference per time step between {case_reference} and {case} is {norm_average}.')
for case, result in results.items():
    iterations = np.array(result['iterations']).mean()
    time = result['run_time']
    tools.print_info(f'{case}: average # iterations = {iterations:0.2f} and elapsed time = {time:0.3f}s')
    tools.print_info('\t', result['iterations'])
