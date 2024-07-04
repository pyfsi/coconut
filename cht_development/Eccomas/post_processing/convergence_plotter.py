from coconut.tools import flatten_concatenation

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import pickle
from scipy import integrate

debug = True

# different cases to be plotted
common_path = '../grid_convergence/'
if debug:
    case_paths = ['M3/case_results.pickle']
    legend_entries = ['M3']
else:
    case_paths = ['M1/case_results.pickle', 'M2/case_results.pickle', 'M3/case_results.pickle', 'M4/case_results.pickle']
    legend_entries = ['M1', 'M2', 'M3', 'M4']

# load cases
results = {}
for name, path in zip(legend_entries, case_paths):
    with open(os.path.join(common_path, path), 'rb') as file:
        results.update({name: pickle.load(file)})

# make figures
line_styles = ['b*', 'r*', 'g*', 'm*', 'y*', 'k*', 'k*', 'k*', 'k*', 'k*', 'k*', 'k*', 'k*']
lines_temperature = []
lines_heat_flux = []
for sol, itf, var, uni in (('solution_x', 'interface_x', 'temperature', 'K'), ('solution_y', 'interface_y', 'heat_flux', 'W/m^2')):
    for j, name in enumerate(legend_entries):
        interface = results[name][itf]
        # only works if interface consists of a single model part!
        model_part = interface.get_model_part(interface.get_model_part_names()[0])
        coordinates = np.zeros((model_part.size, 3))
        for i, direction in enumerate(['x0', 'y0', 'z0']):
            coordinates[:, i] = getattr(model_part, direction)

        if var == 'temperature':
            for k in range(8):
                solution = results[name][sol][:,k]
                line, = plt.plot(coordinates[:,0], solution, line_styles[k], label='ts ' + str(k+1))
                lines_temperature.append(line)
        else:
            for k in range(8):
                solution = results[name][sol][:,k]
                line, = plt.plot(coordinates[:, 0], solution, line_styles[k], label='ts ' + str(k+1))
                lines_heat_flux.append(line)

    # Plot temperature
    if var == 'temperature':
        plt.ylabel('Interface temperature [K]')
        plt.xlabel('Location [m]')
        plt.legend(handles=lines_temperature)
        plt.savefig('itf-temp-M3.svg')
        plt.show()
        plt.close()

    # Plot heat flux
    if var == 'heat_flux':
        plt.ylabel('Interface heat flux [W/m^2]')
        plt.xlabel('Location [m]')
        plt.legend(handles=lines_heat_flux)
        plt.savefig('itf-hf-M3.svg')
        plt.show()
        plt.close()