from coconut import tools

import numpy as np
import os
import matplotlib.pyplot as plt
import pickle

# This script plots the convergence residuals of the coupling iterations in every time step.
# To generate a result file, include a non-zero int 'save_results' in the settings of the coupled solver.
# Give a name to the case by including the string {'case_name': 'a_name'} in the settings of the coupled solver.

tolerance = 1e-8  # tolerance

# different cases to be plotted
common_path = '../../examples/'
case_paths = ['tube_tube_flow_tube_structure/case_results.pickle']
legend_entries = ['case_results']

# load cases
results = {}
for name, path in zip(legend_entries, case_paths):
    with open(os.path.join(common_path, path), 'rb') as file:
        results.update({name: pickle.load(file)})

# reference case
case_reference = legend_entries[0]


def zero_to_nan(values):
    """Replace every 0 with 'nan' and return a copy."""
    return [float('nan') if x == 0 else x for x in values]


# find maximum number of iterations per time step over all cases
it_max = 0
residual_list = []  # list containing residual lists corresponding for the different cases
# each residual list contains a nested list: [ts0[it0, it1, ...], ts1[it0, ...], ...]
for case in legend_entries:
    residual_list.append(results[case]['residual'])
    for ls in residual_list[-1]:
        if len(ls) > it_max:
            it_max = len(ls)

# make figure
plt.figure()
residual = None
for case, residuals in zip(legend_entries, residual_list):
    residual = np.zeros((len(residuals), it_max))
    for i, ls in enumerate(residuals):
        residual[i, :len(ls)] = np.array(ls)
    plt.plot(zero_to_nan(residual.flatten()), 'o-', label=case)
plt.plot(tolerance * np.ones_like(residual.flatten()), 'k')
plt.yscale('log')
plt.ylabel('norm of residual')
plt.xlabel('iteration')
plt.legend()

# average number of iteration per time step
for case, residuals in zip(legend_entries, residual_list):
    iterations = []
    for ls in residuals:
        iterations.append(len(ls))
    avg_iterations = np.array(iterations).mean()
    tools.print_info(f'{case}: average # iterations = {avg_iterations:0.2f}')
    tools.print_info('\t', iterations)

plt.show()
