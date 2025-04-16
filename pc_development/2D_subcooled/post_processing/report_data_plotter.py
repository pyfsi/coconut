import matplotlib.pyplot as plt
import numpy as np
import math as M
import pandas as pd
from scipy.optimize import fsolve
from scipy.special import erf
from scipy.interpolate import interp1d
from coconut.examples.post_processing.post_processing import *

# Start delay of simulation due to starting with LF = 0.01 (based on conduction only Fluent simulation)
t_delay = 12.47 # s

# Read and load coconut data
case_dir = '../Faden_full_5/'
report_file_name = 'report-file.out'
data_solid = np.loadtxt(case_dir + 'CFD_1/' + report_file_name, delimiter=' ', skiprows=3)
data_liquid = np.loadtxt(case_dir + 'CFD_2/' + report_file_name, delimiter=' ', skiprows=3)

# Read and load Fluent comparison data
fluent_dir = 'fluent_reports/'
fluent_report_name = 'report-file.out'
data_fluent = np.loadtxt(fluent_dir + fluent_report_name, delimiter=' ', skiprows=3)

# Construct coconut arrays
area = 0.04 * 1.0 # m^2

# Solid
time = data_solid[:,7] + t_delay
i_end = np.size(time)

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
temp_u1 = np.where(temp_u1_l >= 311.13, np.minimum(temp_u1_s, temp_u1_l), np.maximum(temp_u1_s, temp_u1_l))
temp_u2 = np.where(temp_u2_l >= 311.13, np.minimum(temp_u2_s, temp_u2_l), np.maximum(temp_u2_s, temp_u2_l))
temp_u3 = np.where(temp_u3_l >= 311.13, np.minimum(temp_u3_s, temp_u3_l), np.maximum(temp_u3_s, temp_u3_l))
temp_l1 = np.where(temp_l1_l >= 311.13, np.minimum(temp_l1_s, temp_l1_l), np.maximum(temp_l1_s, temp_l1_l))
temp_l2 = np.where(temp_l2_l >= 311.13, np.minimum(temp_l2_s, temp_l2_l), np.maximum(temp_l2_s, temp_l2_l))

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
line_sim, = plt.plot(time[ts_0:], LF[ts_0:], '-k', label='Coconut')
line_fl_1, = plt.plot(time_fl[ts_0_fl:], LF_1_fl[ts_0_fl:], '--r', label='Fluent 1')
#line_fl_2, = plt.plot(time_fl[ts_0_fl:], LF_2_fl[ts_0_fl:], '--b', label='Fluent 2')
plt.ylabel('Liquid fraction [-]')
plt.xlabel('Time [s]')
plt.legend(handles=[line_sim, line_fl_1])
plt.tight_layout()
plt.savefig('Report_figures/liquid_fraction.png')
plt.show()
plt.close()

# Max. velocity magnitude
line_sim, = plt.plot(time[ts_0:], v_max[ts_0:], '-k', label='Coconut')
line_fl, = plt.plot(time_fl[ts_0_fl:], v_max_fl[ts_0_fl:], '--r', label='Fluent')
plt.ylabel('Max. velocity mag. [m/s]')
plt.xlabel('Time [s]')
plt.legend(handles=[line_sim, line_fl])
plt.tight_layout()
plt.savefig('Report_figures/max_velocity.png')
plt.show()
plt.close()

# Heat flux @ heated wall
line_sim, = plt.plot(time[ts_0:], q_heat[ts_0:], '-k', label='Coconut')
line_fl, = plt.plot(time_fl[ts_0_fl:], q_heat_fl[ts_0_fl:], '--r', label='Fluent')
plt.ylabel('Area avg. heat flux [W/m^2]')
plt.xlabel('Time [s]')
plt.legend(handles=[line_sim, line_fl])
plt.tight_layout()
plt.savefig('Report_figures/hf_heated_wall.png')
plt.show()
plt.close()

# Heat flux @ cooled wall
line_sim, = plt.plot(time[ts_0:], q_cool[ts_0:], '-k', label='Coconut')
line_fl, = plt.plot(time_fl[ts_0_fl:], q_cool_fl[ts_0_fl:], '--r', label='Fluent')
plt.ylabel('Area avg. heat flux [W/m^2]')
plt.xlabel('Time [s]')
plt.legend(handles=[line_sim, line_fl])
plt.tight_layout()
plt.savefig('Report_figures/hf_cooled_wall.png')
plt.show()
plt.close()

# Temperature plot (u1)
line_sim, = plt.plot(time[ts_0:], temp_u1[ts_0:], '-k', label='Coconut')
line_fl, = plt.plot(time_fl[ts_0_fl:], temp_u1_fl[ts_0_fl:], '--r', label='Fluent')
plt.ylabel('Temperature (u1) [K]')
plt.xlabel('Time [s]')
plt.legend(handles=[line_sim, line_fl])
plt.tight_layout()
plt.savefig('Report_figures/Temperature_u1.png')
plt.show()
plt.close()

# Temperature plot (u2)
line_sim, = plt.plot(time[ts_0:], temp_u2[ts_0:], '-k', label='Coconut')
line_fl, = plt.plot(time_fl[ts_0_fl:], temp_u2_fl[ts_0_fl:], '--r', label='Fluent')
plt.ylabel('Temperature (u2) [K]')
plt.xlabel('Time [s]')
plt.legend(handles=[line_sim, line_fl])
plt.tight_layout()
plt.savefig('Report_figures/Temperature_u2.png')
plt.show()
plt.close()

# Temperature plot (u3)
line_sim, = plt.plot(time[ts_0:], temp_u3[ts_0:], '-k', label='Coconut')
line_fl, = plt.plot(time_fl[ts_0_fl:], temp_u3_fl[ts_0_fl:], '--r', label='Fluent')
plt.ylabel('Temperature (u3) [K]')
plt.xlabel('Time [s]')
plt.legend(handles=[line_sim, line_fl])
plt.tight_layout()
plt.savefig('Report_figures/Temperature_u3.png')
plt.show()
plt.close()

# Temperature plot (l1)
line_sim, = plt.plot(time[ts_0:], temp_l1[ts_0:], '-k', label='Coconut')
line_fl, = plt.plot(time_fl[ts_0_fl:], temp_l1_fl[ts_0_fl:], '--r', label='Fluent')
plt.ylabel('Temperature (l1) [K]')
plt.xlabel('Time [s]')
plt.legend(handles=[line_sim, line_fl])
plt.tight_layout()
plt.savefig('Report_figures/Temperature_l1.png')
plt.show()
plt.close()

# Temperature plot (l2)
line_sim, = plt.plot(time[ts_0:], temp_l2[ts_0:], '-k', label='Coconut')
line_fl, = plt.plot(time_fl[ts_0_fl:], temp_l2_fl[ts_0_fl:], '--r', label='Fluent')
plt.ylabel('Temperature (l2) [K]')
plt.xlabel('Time [s]')
plt.legend(handles=[line_sim, line_fl])
plt.tight_layout()
plt.savefig('Report_figures/Temperature_l2.png')
plt.show()
plt.close()