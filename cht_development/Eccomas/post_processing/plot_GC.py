import numpy as np
import os
import matplotlib.pyplot as plt
import pickle

debug = False

# different cases to be plotted
common_path = '../grid_convergence/'
if debug:
    case_paths = ['M1/case_results.pickle', 'M2/case_results.pickle']
    legend_entries = ['M1', 'M2']
else:
    case_paths = ['M1/case_results.pickle', 'M2/case_results.pickle', 'M3/case_results.pickle', 'M4/case_results.pickle']
    legend_entries = ['M1', 'M2', 'M3', 'M4']

# load cases
results = {}
for name, path in zip(legend_entries, case_paths):
    with open(os.path.join(common_path, path), 'rb') as file:
        results.update({name: pickle.load(file)})

for name in legend_entries:
    print(results[name]['info'])
    print('Run time:', results[name]['run_time']/60 ,'minutes')

# make figures
line_styles = ['b--', 'r--', 'g--', 'm--']
lines_temperature = []
lines_heat_flux = []
for sol, itf, var, uni in (('solution_x', 'interface_x', 'temperature', 'K'), ('solution_y', 'interface_y', 'heat_flux', 'W/m^2')):
    for j, name in enumerate(legend_entries):
        interface = results[name][itf]
        # only works if interface consists of a single model part!
        model_part = interface.get_model_part(interface.get_model_part_names()[0])
        x0 = getattr(model_part, 'x0')
        args = np.unique(x0, return_index=True)[1].tolist()
        x = x0[args]

        if var == 'temperature':
            solution = results[name][sol][args,-1]
            line, = plt.plot(x, solution, line_styles[j], label=name)
            lines_temperature.append(line)
        else:
            solution = results[name][sol][args,-1]
            line, = plt.plot(x, solution, line_styles[j], label=name)
            lines_heat_flux.append(line)

    # Plot temperature
    if var == 'temperature':
        plt.ylabel('Interface temperature [K]')
        plt.xlabel('Location [m]')
        plt.legend(handles=lines_temperature)
        plt.savefig('./figures/itf-temp-GC.svg')
        plt.show()
        plt.close()

    # Plot heat flux
    if var == 'heat_flux':
        plt.ylabel('Interface heat flux [W/m^2]')
        plt.xlabel('Location [m]')
        plt.legend(handles=lines_heat_flux)
        plt.savefig('./figures/itf-hf-GC.svg')
        plt.show()
        plt.close()