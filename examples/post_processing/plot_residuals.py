from coconut import tools

import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
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

# find maximum number of iterations per time step over all cases
it_max = 0
residual_list = []  # list containing residual lists corresponding for the different cases
iterations_list = []  # list containing number of iterations corresponding for the different cases
# each residual list contains a nested list: [ts0[it0, it1, ...], ts1[it0, ...], ...]
for case in legend_entries:
    residual_list.append(results[case]['residual'])
    iterations_list.append(results[case]['iterations'])
    for ls in residual_list[-1]:
        if len(ls) > it_max:
            it_max = len(ls) + 1

# make figure
plt.figure()
residual = None
for case, residuals in zip(legend_entries, residual_list):
    residual = np.ones((len(residuals), it_max)) * np.nan
    for i, ls in enumerate(residuals):
        residual[i, :len(ls)] = np.array(ls)
    plt.plot(residual.flatten(), 'o-', label=case)
plt.axhline(tolerance, color='k', label=f'tolerance {tolerance}')
plt.gca().xaxis.set_major_formatter(FuncFormatter(lambda value, _: int(value // it_max + 1)))
plt.yscale('log')
plt.ylabel('norm of residual')
plt.xlabel('time step')
plt.legend()

# average number of iteration per time step
for case, iterations in zip(legend_entries, iterations_list):
    avg_iterations = np.mean(iterations)
    tools.print_info(f'{case}: average # iterations = {avg_iterations:0.2f}')
    tools.print_info('\t', iterations)

plt.show()
