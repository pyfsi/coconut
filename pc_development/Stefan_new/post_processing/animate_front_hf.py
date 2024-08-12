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

itf_faces = 10

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
            for j in [0 + i*3 for i in range((np.shape(solution)[0]-itf_faces)//3)]:
                disp_x += solution[j,:]
            disp_x /= (np.shape(solution)[0]-itf_faces)//3
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

# Plot liquid fraction
L = 0.05 # m, length of original solid domain
LF = 0.5-1*disp_x/L

# Read Fluent validation files
read_file_1 = pd.read_csv(r'lf_time.out', delimiter='\s+', skiprows=[0, 1, 2]) # Path of Fluent out-file
read_file_1.to_csv(r'lf_time.csv', index=None)
data_array_1 = np.loadtxt('lf_time.csv', delimiter=',')
LF_val_1 = data_array_1[:,1]
time_val_1 = data_array_1[:,2]
try:
    os.remove("lf_time.csv")
except:
    pass

# Simple steady model
k = 1.5 # W/mK, conduction coefficient of liquid
dT = 10 # K, the maintained temperature difference over the liquid domain
rho = 870 # kg/m^3, PCM density
LH = 179000 # J/kg, latent heat
B = (k*dT)/(rho*LH) # m^2/s

time_ana = np.linspace(0, 3600, 3601)
dx = -0.5*L + 0.5*np.sqrt(L**2 + 4*B*time_val_1)
LF_ana = 0.5 + dx/L

line1, = plt.plot(time, LF.flatten(), '-k', label="CoCoNuT")
line2, = plt.plot(time_val_1, LF_val_1, '--r', label="Fluent")
line3, = plt.plot(time_val_1, LF_ana, '--b', label="Ana - steady")
plt.ylabel('Liquid fraction')
plt.xlabel('Time [s]')
plt.legend(handles=[line1, line2, line3])
plt.savefig('liq-frac.png')
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