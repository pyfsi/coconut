import matplotlib.pyplot as plt
import numpy as np
import math as M
import pandas as pd
from coconut.examples.post_processing.post_processing import *

# ----- Settings -----
temp_plots = False
faden = True
faden_dir = 'Faden_paper/'

# Start delay for LF
t_delay_001, t_delay_002 = 12.47, 61.23
start_lf_001 = False
t_delay = t_delay_001 if start_lf_001 else t_delay_002

# Case directories and report files
case_dir = '../Faden_split_4/'
report_file_name = 'report-file.out'
fluent_dir = 'fluent_reports/'
fluent_report_name = 'report-file-temp.out'

# ----- Load data -----
data_solid = np.loadtxt(f'{case_dir}CFD_1/{report_file_name}', delimiter=' ', skiprows=3)
data_liquid = np.loadtxt(f'{case_dir}CFD_2/{report_file_name}', delimiter=' ', skiprows=3)
data_fluent = np.loadtxt(f'{fluent_dir}{fluent_report_name}', delimiter=' ', skiprows=3)

# ----- Construct arrays -----
area = 0.04 * 0.08
area_2D = 0.04 * 1
full_vol = 0.04 * 0.04

time = data_solid[:,7] + t_delay
time_l = data_liquid[:,9] + t_delay
i_end = min(len(time), len(time_l))

# Solid
q_cool = -area * data_solid[:i_end,1]
temp_s = data_solid[:i_end,2:7]  # u1,u2,u3,l1,l2

# Liquid
vol_liquid = data_liquid[:i_end,1]
v_max = data_liquid[:i_end,2]
q_heat = area * data_liquid[:i_end,3]
temp_l = data_liquid[:i_end,4:9]

# Liquid fraction
LF = vol_liquid / full_vol

# Merge solid/liquid temperatures
tol = 3
temp_combined = np.where(temp_l >= 311.13 - tol, np.minimum(temp_s, temp_l), np.maximum(temp_s, temp_l))
temp_u1, temp_u2, temp_u3, temp_l1, temp_l2 = temp_combined.T

# Fluent arrays
solid_mass_1, solid_mass_2 = data_fluent[:,1], data_fluent[:,2]
v_max_fl = data_fluent[:,3]
q_heat_fl = area * data_fluent[:,4]
q_cool_fl = -area * data_fluent[:,5]
temp_fl = data_fluent[:,6:11]  # u1,u2,u3,l1,l2
time_fl = data_fluent[:,11]

mass_ini = 1.3886624 # kg
LF_1_fl = 1 - (solid_mass_1 / mass_ini)
LF_2_fl = 1 - (solid_mass_2 / mass_ini)

# ----- Energy balance -----
net_heat_flux = area_2D * (data_fluent[:,4] + data_fluent[:,5])
heat_flux_int = np.trapz(net_heat_flux, x=time_fl)

mass_end = 1.3157895 # kg
h_hot = 265038.54 # J/kg (enthalpy at T_m)
H_end = 149636.7 # J
H_ini = -18108.53 # J

# INTERMEDIATE VALUES OF NEW SIM!! ALSO SEE REPORT FILE OF FLUENT!!
mass_end = 1.3527841 # kg
h_hot = 265038.54 # J/kg (enthalpy at T_m)
H_end = 70626.203 # J
H_ini = -18108.53 # J

dH = H_end - H_ini
H_out = (mass_ini - mass_end) * h_hot
heat_net_in = heat_flux_int - H_out
deficit = abs(heat_net_in - dH)

print(f"Total heat transfer at walls: {heat_flux_int:.3f} J")
print(f"Total energy lost at outlet: {H_out:.3f} J")
print(f"Net added energy: {heat_net_in:.3f} J")
print(f"Enthalpy increase: {dH:.3f} J")
print(f"Energy deficit: {(deficit / dH) * 100:.3f} %")

# ----- Helper function for plots -----
def plot_with_faden(time_sim, y_sim, label_sim, time_fl, y_fl, label_fl,
                    faden_files=None, faden_labels=None, line_styles=None,
                    ylabel='', xlabel='Time [s]', save_name='plot.png', scale_time=False,
                    exp_range_files=None, y_min=None, y_max=None):
    """
    Plots simulation, Fluent, and Faden data with experimental ranges as discrete error bars.

    Parameters:
    - exp_range_files: list of tuples (file, label) for experimental range CSVs
    - y_min, y_max: optional floats to set y-axis limits
    """
    # Partitioned simulation
    plt.plot(time_sim, y_sim, '-k', label=label_sim)

    # Fluent
    plt.plot(time_fl, y_fl, '--r', label=label_fl)

    # Faden numerical/experimental points (not ranges)
    if faden and faden_files:
        for j, file in enumerate(faden_files):
            t_f, y_f = np.loadtxt(f'{faden_dir}{file}', skiprows=1, delimiter=',', unpack=True)
            if scale_time:
                t_f *= 60
            plt.plot(t_f, y_f, line_styles[j], label=faden_labels[j])

    # Faden experimental ranges as discrete error bars
    if faden and exp_range_files:
        for file, label in exp_range_files:
            data = np.loadtxt(f'{faden_dir}{file}', delimiter=',', skiprows=1)
            time_exp, LF_exp = data[:, 0], data[:, 1]
            unique_times = np.unique(time_exp)
            LF_min, LF_max = [], []

            for t in unique_times:
                LF_at_t = LF_exp[time_exp == t]
                LF_min.append(np.min(LF_at_t))
                LF_max.append(np.max(LF_at_t))

            LF_min = np.array(LF_min)
            LF_max = np.array(LF_max)

            # Compute mean and error
            LF_mean = (LF_max + LF_min) / 2
            LF_err = np.vstack((LF_mean - LF_min, LF_max - LF_mean))

            # Small offset if min==max for visibility
            LF_err[LF_err == 0] = 0.002

            if scale_time:
                unique_times = unique_times * 60

            plt.errorbar(unique_times, LF_mean, yerr=LF_err, fmt='g.', capsize=3, label=label)

    plt.xlabel(xlabel, fontsize=16)
    plt.ylabel(ylabel, fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    if y_min is not None or y_max is not None:
        plt.ylim(y_min, y_max)
    plt.legend(fontsize=14)
    plt.tight_layout()
    plt.savefig(f'Report_figures/{save_name}', dpi=150)
    plt.show()
    plt.close()


# ----- Plot main results for paper -----
ts_0, ts_0_fl = M.ceil(1.0/0.1), 118

# Liquid fraction
plot_with_faden(
    time[ts_0:], LF[ts_0:], 'Partitioned',
    time_fl[ts_0_fl:], LF_1_fl[ts_0_fl:], 'Fixed grid',
    faden_files=['Faden_LF_sim.csv'],                  # numerical points
    faden_labels=['Faden - num.'],
    line_styles=['b.'],
    ylabel='Liquid fraction [-]',
    save_name='liquid_fraction.png',
    scale_time=True,
    exp_range_files=[('Faden_LF_exp_itv.csv', 'Faden - exp.')])


# Heat flux - heated wall (with Faden data)
plot_with_faden(time[ts_0:], q_heat[ts_0:], 'Partitioned',
                time_fl[ts_0_fl:], q_heat_fl[ts_0_fl:], 'Fixed grid',
                faden_files=['Faden_HF_heated_exp.csv', 'Faden_HF_heated_sim.csv'],
                faden_labels=['Faden - exp.', 'Faden - num.'],
                line_styles=['g--', 'b--'],
                ylabel='Heat transfer rate [W]',
                save_name='hf_heated_wall.png',
                y_min=2, y_max=4)

# Heat flux - cooled wall (with Faden data)
plot_with_faden(time[ts_0:], q_cool[ts_0:], 'Partitioned',
                time_fl[ts_0_fl:], q_cool_fl[ts_0_fl:], 'Fixed grid',
                faden_files=['Faden_HF_cooled_exp.csv', 'Faden_HF_cooled_sim.csv'],
                faden_labels=['Faden - exp.', 'Faden - num.'],
                line_styles=['g--', 'b--'],
                ylabel='Heat transfer rate [W]', save_name='hf_cooled_wall.png')


# ----- Other plots remain intact -----
# Max. velocity
plt.plot(time[ts_0:], v_max[ts_0:]*1000, '-k', label='Partitioned')
plt.plot(time_fl[ts_0_fl:], v_max_fl[ts_0_fl:]*1000, '--r', label='Fixed grid')
plt.ylabel('Max. velocity mag. [mm/s]', fontsize=16)
plt.xlabel('Time [s]', fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(fontsize=16)
plt.tight_layout()
plt.savefig('Report_figures/max_velocity.png', dpi=150)
plt.show()
plt.close()

# Temperature plots
if temp_plots:
    labels = ['u1','u2','u3','l1','l2']
    for i, lbl in enumerate(labels):
        plt.plot(time[ts_0:], temp_combined[ts_0:,i], '-k', label='Coconut')
        plt.plot(time_fl[ts_0_fl:], temp_fl[ts_0_fl:,i], '--r', label='Fluent')
        plt.ylabel(f'Temperature ({lbl}) [K]')
        plt.xlabel('Time [s]')
        plt.legend()
        plt.tight_layout()
        plt.savefig(f'Report_figures/Temperature_{lbl}.png')
        plt.show()
        plt.close()
