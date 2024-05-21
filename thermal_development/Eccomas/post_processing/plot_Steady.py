import numpy as np
import os
import matplotlib.pyplot as plt
import pickle

debug = False
v = 36

# different cases to be plotted
dict = {2.4: '2-4m_s', 12: '12m_s', 24: '24m_s', 36: '36m_s'}
common_path = '../../../thermal_development/Eccomas/Steady/' + dict[2.4] + '/'
f = '/case_results.pickle'

if debug:
    case_paths = ['FFTB_relaxation' + f]
    legend_entries = ['FFTB_relaxation']
else:
    case_paths = ['TFFB_relaxation' + f, 'FFTB_relaxation' + f, 'TFFB_aitken' + f, 'FFTB_aitken' + f, 'TFFB_iqni' + f, 'FFTB_iqni' + f]
    legend_entries = ['TFFB_relaxation', 'FFTB_relaxation', 'TFFB_aitken', 'FFTB_aitken', 'TFFB_iqni', 'FFTB_iqni']

# load cases
results = {}
for name, path in zip(legend_entries, case_paths):
    with open(os.path.join(common_path, path), 'rb') as file:
        results.update({name: pickle.load(file)})

for name in legend_entries:
    print(results[name]['info'])
    print('Run time:', results[name]['run_time']/60 ,'minutes')

# make figures
line_styles = ['b--', 'b.-', 'g--', 'g.-', 'r--', 'r.-']
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
            if np.mean(solution) < 0:
                solution = -1*solution
            line, = plt.plot(x, solution, line_styles[j], label=name)
            lines_heat_flux.append(line)

    # Plot temperature
    if var == 'temperature':
        plt.ylabel('Interface temperature [K]')
        plt.xlabel('Location [m]')
        plt.legend(handles=lines_temperature)
        plt.title(dict[v])
        plt.savefig('./figures/itf-temp-steady-' + dict[v] + '.svg')
        plt.show()
        plt.close()

    # Plot heat flux
    if var == 'heat_flux':
        plt.ylabel('Interface heat flux [W/m^2]')
        plt.xlabel('Location [m]')
        plt.legend(handles=lines_heat_flux)
        plt.title(dict[v])
        plt.savefig('./figures/itf-hf-steady-' + dict[v] + '.svg')
        plt.show()
        plt.close()