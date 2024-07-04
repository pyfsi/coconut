from coconut.tools import flatten_concatenation

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import pickle
from scipy import integrate

# different cases to be plotted
common_path = '../'
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
            avg_T = np.array([np.mean(solution[:,i]) for i in range(np.shape(solution)[1])])
            # avg_T[0] = 283.15 # set manually
            time = time_step_start + dt*np.array(range(np.shape(solution)[1]))

        if var == "heat_flux":
            avg_hf = np.array([np.mean(solution[:, i]) for i in range(np.shape(solution)[1])])

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
plt.savefig('itf-temp-forced.svg')
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

# Energy balance
read_file_sol = pd.read_csv(os.path.join(common_path, 'CFD_1/report-file.out'), delimiter='\s+', skiprows=[0, 1, 2])
read_file_liq = pd.read_csv(os.path.join(common_path, 'CFD_2/report-file.out'), delimiter='\s+', skiprows=[0, 1, 2])
read_file_sol.to_csv(os.path.join(common_path, 'CFD_1/report-file.csv'), index=None)
read_file_liq.to_csv(os.path.join(common_path, 'CFD_2/report-file.csv'), index=None)

data_sol = np.loadtxt(os.path.join(common_path, 'CFD_1/report-file.csv'), delimiter=',')
data_liq = np.loadtxt(os.path.join(common_path, 'CFD_2/report-file.csv'), delimiter=',')

# Exclude 1st second
data_sol = data_sol[10:,:]
data_liq = data_liq[10:,:]

temp_sol = data_sol[:,1]
itf_flux_sol = data_sol[:,2]
wall_flux_sol = data_sol[:,3]
time_s = data_sol[:,4]

temp_liq = data_liq[:,1]
itf_flux_liq = data_liq[:,2]
wall_flux_liq = data_liq[:,3]
integral_out = data_liq[:,4]
time_l = data_liq[:,5]

line5, = plt.plot(time_s, -1*itf_flux_sol, 'b', label='Solid')
line6, = plt.plot(time_l, itf_flux_liq, 'g', label='Liquid')
line7, = plt.plot(time[10:], -1*avg_hf[10:], 'r', label='Coconut')
plt.ylabel('Heat flux [W/m^2]')
plt.xlabel('Time [s]')
plt.legend(handles=[line5, line6, line7])
plt.show()
plt.close()

# Properties and dimensions
rho_l = 780 # kg/m^3
rho_s = 870 # kg/m^3
cp = 2500 # J/kgK
V = 0.05*0.05*1 # m^3
A = 0.05*1 # m^2
T_ini_sol = 273.15 + 30 # °C
T_ini_liq = 273.15 + 10 # °C
v = 0.005 # m/s

# Calculation solid
htr_wall_sol_int = integrate.trapezoid(wall_flux_sol*A, time_s)
htr_itf_sol_int = integrate.trapezoid(itf_flux_sol*A, time_s)

SH_sol = rho_s*V*cp*(temp_sol[-1]-temp_sol[0])

sol_balance = SH_sol/(htr_itf_sol_int+htr_wall_sol_int)
print('Solid balance is', sol_balance*100, '%.')

# Calculation liquid
delta_t = time_l[1:] - time_l[:-1]
delta_T = temp_liq[1:] - temp_liq[:-1]
dT_dt = np.divide(delta_T, delta_t)
net_flow = rho_l*A*v*cp*T_ini_liq - (integral_out[1:] + integral_out[:-1])/2
flux_liq_avg = ((itf_flux_liq+wall_flux_liq)[1:] + (itf_flux_liq+wall_flux_liq)[:-1])/2

Q_dot = A*flux_liq_avg
delta_SH_liq = rho_l*V*cp*dT_dt

line8, = plt.plot(time_l[:-1], Q_dot, 'b', label='Q_dot')
line9, = plt.plot(time_l[:-1], delta_SH_liq - net_flow, 'g', label='SH rate + net outflow')
plt.ylabel('Heat rate [W]')
plt.xlabel('Time [s]')
plt.legend(handles=[line8, line9])
plt.show()
plt.close()

flux_tot_liq = integrate.trapezoid((itf_flux_liq+wall_flux_liq)*A, time_l)
htr_wall_liq = integrate.trapezoid(wall_flux_liq*A, time_l)
outflow_tot_liq = integrate.trapezoid(rho_l*A*v*cp*T_ini_liq - integral_out, time_l)
SH_liq = rho_l*V*cp*(temp_liq[-1] - temp_liq[0])
liq_balance = (SH_liq-outflow_tot_liq)/flux_tot_liq
print('Liquid balance is', liq_balance*100, '%.')

global_balance = (SH_sol+SH_liq)/(htr_wall_liq+htr_wall_sol_int+outflow_tot_liq)
print('Global balance is', global_balance*100, '%.')

# Cleanup
os.remove(os.path.join(common_path, 'CFD_1/report-file.csv'))
os.remove(os.path.join(common_path, 'CFD_2/report-file.csv'))