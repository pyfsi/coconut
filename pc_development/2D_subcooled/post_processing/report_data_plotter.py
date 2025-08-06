import matplotlib.pyplot as plt
import numpy as np
import math as M
import pandas as pd
from coconut.examples.post_processing.post_processing import *

# Plot temperatures?
temp_plots = False

# Plot Faden validation data?
faden = True
faden_dir = 'Faden_paper/'

# Start delay of simulation due to starting with LF = 0.01 & 0.02 (based on conduction only Fluent simulation)
t_delay_001 = 12.47 # s
t_delay_002 = 61.23 # s 61.23 or 58.79

# Read and load coconut data
start_lf_001 = False
case_dir = '../Faden_split_3_coarse/'
report_file_name = 'report-file.out'
data_solid = np.loadtxt(case_dir + 'CFD_1/' + report_file_name, delimiter=' ', skiprows=3)
data_liquid = np.loadtxt(case_dir + 'CFD_2/' + report_file_name, delimiter=' ', skiprows=3)

# Read and load Fluent comparison data
fluent_dir = 'fluent_reports/'
fluent_report_name = 'report-file-new.out'
data_fluent = np.loadtxt(fluent_dir + fluent_report_name, delimiter=' ', skiprows=3)

t_delay = t_delay_001 if start_lf_001 else t_delay_002

# Construct coconut arrays
area = 0.04 * 1.0 # m^2
time = data_solid[:,7]
time_l = data_liquid[:,9]
i_end = min(np.size(time), np.size(time_l))

# Solid
time = data_solid[:i_end,7] + t_delay
q_cool = -1 * area * data_solid[:i_end,1] # [W]
temp_u1_s = data_solid[:i_end,2]
temp_u2_s = data_solid[:i_end,3]
temp_u3_s = data_solid[:i_end,4]
temp_l1_s = data_solid[:i_end,5]
temp_l2_s = data_solid[:i_end,6]

# Liquid
vol_liquid = data_liquid[:i_end,1]
v_max = data_liquid[:i_end,2]
q_heat = area * data_liquid[:i_end,3] # [W]
temp_u1_l = data_liquid[:i_end,4]
temp_u2_l = data_liquid[:i_end,5]
temp_u3_l = data_liquid[:i_end,6]
temp_l1_l = data_liquid[:i_end,7]
temp_l2_l = data_liquid[:i_end,8]

# Calculate LF coconut
full_vol = 0.04 * 0.04 # m^3
LF = vol_liquid / full_vol

# Merge temperature arrays
tol = 2
temp_u1 = np.where(temp_u1_l >= 311.13 - tol, np.minimum(temp_u1_s, temp_u1_l), np.maximum(temp_u1_s, temp_u1_l))
temp_u2 = np.where(temp_u2_l >= 311.13 - tol, np.minimum(temp_u2_s, temp_u2_l), np.maximum(temp_u2_s, temp_u2_l))
temp_u3 = np.where(temp_u3_l >= 311.13 - tol, np.minimum(temp_u3_s, temp_u3_l), np.maximum(temp_u3_s, temp_u3_l))
temp_l1 = np.where(temp_l1_l >= 311.13 - tol, np.minimum(temp_l1_s, temp_l1_l), np.maximum(temp_l1_s, temp_l1_l))
temp_l2 = np.where(temp_l2_l >= 311.13 - tol, np.minimum(temp_l2_s, temp_l2_l), np.maximum(temp_l2_s, temp_l2_l))

# Construct fluent arrays
solid_mass_1 = data_fluent[:,1]
solid_mass_2 = data_fluent[:,2]
v_max_fl = data_fluent[:,3]
q_heat_fl = area * data_fluent[:,4] # [W]
q_cool_fl = -1 * area * data_fluent[:,5] # [W]
temp_u1_fl = data_fluent[:,6]
temp_u2_fl = data_fluent[:,7]
temp_u3_fl = data_fluent[:,8]
temp_l1_fl = data_fluent[:,9]
temp_l2_fl = data_fluent[:,10]
time_fl = data_fluent[:,11]

# Calculate LF Fluent
mass_ini = 1.388662402342756 # kg
LF_1_fl = 1 - (solid_mass_1 / mass_ini)
LF_2_fl = 1 - (solid_mass_2 / mass_ini)

# Plots (simulation data starting from t_0)
ts_0 = M.ceil(1.0/0.1)
ts_0_fl = 118

# Create plots

# Liquid fraction
line_sim, = plt.plot(time[ts_0:], LF[ts_0:], '-k', label='Partitioned')
line_fl_1, = plt.plot(time_fl[ts_0_fl:], LF_1_fl[ts_0_fl:], '--r', label='Fixed grid')
lines_lf = [line_sim, line_fl_1]

# Read and load Faden paper comparison data for liquid fraction
if faden:
    faden_data = ['Faden_LF_exp_itv.csv', 'Faden_LF_sim.csv']
    faden_legend = ['Faden - exp. 0.95 itv', 'Faden - numerical']
    line_styles = ['g.', 'b.']

    for j, file in enumerate(faden_data):
        time_Fa, LF_Fa = np.loadtxt(faden_dir + file, skiprows=1, delimiter=',', unpack=True)
        line, = plt.plot(time_Fa*60, LF_Fa, line_styles[j], label=faden_legend[j]) # time is in minutes
        lines_lf.append(line)

plt.ylabel('Liquid fraction [-]', fontsize=16)
plt.xlabel('Time [s]', fontsize=16)
#plt.xlim((-100, 2500))
#plt.ylim((0, 0.20))
plt.xticks(fontsize=14)  # Increase font size of x-ticks
plt.yticks(fontsize=14)  # Increase font size of y-ticks
plt.legend(handles=lines_lf, fontsize=14)
plt.tight_layout()
plt.savefig('Report_figures/liquid_fraction.png', dpi=150)
plt.show()
plt.close()

# Max. velocity magnitude
line_sim, = plt.plot(time[ts_0:], v_max[ts_0:]*1000, '-k', label='Partitioned')
line_fl, = plt.plot(time_fl[ts_0_fl:], v_max_fl[ts_0_fl:]*1000, '--r', label='Fixed grid')
plt.ylabel('Max. velocity mag. [mm/s]', fontsize=16)
plt.xlabel('Time [s]', fontsize=16)
#plt.xlim((-100, 2500))
plt.xticks(fontsize=14)  # Increase font size of x-ticks
plt.yticks(fontsize=14)  # Increase font size of y-ticks
plt.legend(handles=[line_sim, line_fl], fontsize=16)
plt.tight_layout()
plt.savefig('Report_figures/max_velocity.png', dpi=150)
plt.show()
plt.close()

# Heat flux @ heated wall
line_sim, = plt.plot(time[ts_0:], q_heat[ts_0:], '-k', label='Partitioned')
line_fl, = plt.plot(time_fl[ts_0_fl:], q_heat_fl[ts_0_fl:], '--r', label='Fixed grid')
plt.ylabel('Area avg. heat flux [W/m^2]', fontsize=16)
plt.xlabel('Time [s]', fontsize=16)
#plt.xlim((-100, 2500))
#plt.ylim((0, 120))
plt.xticks(fontsize=14)  # Increase font size of x-ticks
plt.yticks(fontsize=14)  # Increase font size of y-ticks
plt.legend(handles=[line_sim, line_fl], fontsize=16)
plt.tight_layout()
plt.savefig('Report_figures/hf_heated_wall.png', dpi=150)
plt.show()
plt.close()

# Heat flux @ cooled wall
line_sim, = plt.plot(time[ts_0:], q_cool[ts_0:], '-k', label='Partitioned')
line_fl, = plt.plot(time_fl[ts_0_fl:], q_cool_fl[ts_0_fl:], '--r', label='Fixed grid')
plt.ylabel('Area avg. heat flux [W/m^2]', fontsize=16)
plt.xlabel('Time [s]', fontsize=16)
#plt.xlim((-100, 2500))
#plt.ylim((-0.3, 3.5))
plt.xticks(fontsize=14)  # Increase font size of x-ticks
plt.yticks(fontsize=14)  # Increase font size of y-ticks
plt.legend(handles=[line_sim, line_fl], fontsize=16)
plt.tight_layout()
plt.savefig('Report_figures/hf_cooled_wall.png', dpi=150)
plt.show()
plt.close()

# Temperature plots
if temp_plots:
    # Combine temperature data into a single array
    temp_data_sim = np.array([temp_u1, temp_u2, temp_u3, temp_l1, temp_l2]).T  # Shape: (n, 5)
    temp_data_fl = np.array([temp_u1_fl, temp_u2_fl, temp_u3_fl, temp_l1_fl, temp_l2_fl]).T # Shape (n,5)

    # Labels for the plots
    labels = ['u1', 'u2', 'u3', 'l1', 'l2']

    # Loop through each temperature dataset
    for i in range(len(labels)):
        plt.figure() #Creates a new figure for each plot.
        line_sim, = plt.plot(time[ts_0:], temp_data_sim[ts_0:, i], '-k', label='Coconut')
        line_fl, = plt.plot(time_fl[ts_0_fl:], temp_data_fl[ts_0_fl:,i], '--r', label='Fluent')
        plt.ylabel(f'Temperature ({labels[i]}) [K]')
        plt.xlabel('Time [s]')
        plt.legend(handles=[line_sim, line_fl])
        plt.tight_layout()
        plt.savefig(f'Report_figures/Temperature_{labels[i]}.png')
        plt.show()
        plt.close() #Closes the current figure