from coconut.tools import flatten_concatenation

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import pickle

flag = False

# Initial interface temperature
T_ini = 20 # °C

# different cases to be plotted
common_path = '../../../thermal_development/cond-nat_conv/'
case_paths = ['case_results_GS.pickle']
legend_entries = ['case']

# load cases
results = {}
for name, path in zip(legend_entries, case_paths):
    with open(os.path.join(common_path, path), 'rb') as file:
        results.update({name: pickle.load(file)})


print(results['case']['info'])
print('Run time:', results['case']['run_time']/60 ,'minutes')

# make figures
colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']
line_styles = ['-', '--', ':', '-.']
for sol, itf, var, uni in (('solution_x', 'interface_x', 'temperature', 'K'), ('solution_y', 'interface_y', 'heat_flux', 'W/m^2')):
    for j, name in enumerate(legend_entries):
        if var == 'temperature':
            solution = results[name][sol] - 273.15
        else:
            solution = results[name][sol]
        interface = results[name][itf]
        dt = results[name]['delta_t']
        time_step_start = results[name]['timestep_start']

        # Store interface data for temperature plot
        if var == "temperature":
            flag = True
            avg_T = np.array([np.mean(solution[:,i]) for i in range(np.shape(solution)[1])])
            avg_T[0] = T_ini
            time = time_step_start + dt*np.array(range(np.shape(solution)[1]))

if flag:
    # Plot avg interface temperature in time
    # First, read Fluent validation files
    # domain 1: solid
    read_file_1 = pd.read_csv(r'itf-temp-1.out', delimiter='\s+', skiprows=[0, 1, 2]) # Path of Fluent out-file
    read_file_1.to_csv(r'fluent_validation_solid.csv', index=None)
    data_array_1 = np.loadtxt('fluent_validation_solid.csv', delimiter=',')
    T_val_1 = data_array_1[:,1] - 273.15
    time_val_1 = data_array_1[:,2]
    try:
        os.remove("fluent_validation_solid.csv")
    except:
        pass

    # domain 2: liquid
    read_file_2 = pd.read_csv(r'itf-temp-2.out', delimiter='\s+', skiprows=[0, 1, 2])  # Path of Fluent out-file
    read_file_2.to_csv(r'fluent_validation_liquid.csv', index=None)
    data_array_2 = np.loadtxt('fluent_validation_liquid.csv', delimiter=',')
    T_val_2 = data_array_2[:, 1] - 273.15
    time_val_2 = data_array_2[:, 2]
    try:
        os.remove("fluent_validation_liquid.csv")
    except:
        pass

    # Plot
    line1, = plt.plot(time, avg_T, '--r', label="CoCoNuT")
    line2, = plt.plot(time_val_1, T_val_1, '--g', label="Fluent solid")
    line3, = plt.plot(time_val_2, T_val_2, '--b', label="Fluent liquid")
    plt.ylabel('Interface temperature [°C]')
    plt.xlabel('Time [s]')
    plt.ylim((np.min(avg_T)-1, np.max(avg_T)+1))
    plt.legend(handles=[line1, line2, line3])
    plt.savefig('itf-temp-nat.svg')
    plt.show()
    plt.close()

# Plot residuals
res = np.array(flatten_concatenation(results['case']['residual']))
it = np.arange(res.size)
line4, = plt.plot(it, res)
plt.ylabel('Residual')
plt.xlabel('Nr. of iterations')
plt.show()
plt.close()