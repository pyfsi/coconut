from coconut.tools import flatten_concatenation

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import pickle
from scipy import integrate

# different cases to be plotted
common_path = '../../../pc_development/Stefan_new/'
case_paths = ['case_results.pickle']
legend_entries = ['case']

# load cases
results = {}
for name, path in zip(legend_entries, case_paths):
    with open(os.path.join(common_path, path), 'rb') as file:
        results.update({name: pickle.load(file)})


# print(results['case'])

# make figures
colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']
line_styles = ['-', '--', ':', '-.']
for sol, itf, var, uni in (('solution_x', 'interface_x', 'displacement', 'm'), ('solution_y', 'interface_y', 'heat_flux', 'W/m^2')):
    for j, name in enumerate(legend_entries):
        solution = results[name][sol]
        interface = results[name][itf]
        dt = results[name]['delta_t']
        time_step_start = results[name]['timestep_start']

        # Store interface data for temperature plot
        if var == "displacement":
            disp_x = np.zeros((1, np.shape(solution)[1]))
            for j in [0 + i*3 for i in range((np.shape(solution)[0]-1)//3)]:
                disp_x += solution[j,:]
            disp_x /= (np.shape(solution)[0]-1)//3
            #print(disp_x.flatten())

            time = time_step_start + dt*np.array(range(np.shape(solution)[1]))
            #print(time)

        if var == "heat_flux":
            heat_flux = solution[-1,:]
            #print(heat_flux)

# Plot interface displacement en heat flux in time
plt.plot(time, disp_x.flatten(), '-k')
plt.ylabel('Interface displacement (x) [m]')
plt.xlabel('Time [s]')
plt.savefig('itf-disp-x.png')
plt.show()
plt.close()

plt.plot(time, heat_flux, '-k')
plt.ylabel('Heat flux [W/m^2]')
plt.xlabel('Time [s]')
plt.savefig('itf-hf.png')
plt.show()
plt.close()

# Plot residuals
"""
res = np.array(flatten_concatenation(results['case']['residual']))
it = np.arange(res.size)
line4, = plt.plot(it, res)
plt.ylabel('Residual')
plt.xlabel('Nr. of iterations')
plt.show()
plt.close()
"""