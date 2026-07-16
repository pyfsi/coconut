import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
import math
from coconut.examples.post_processing.post_processing import PostProcess

# =============================================================================
# 1. SETTINGS & CONFIGURATION
# =============================================================================
# Directories
SIM_DIR = '../fluent_fluent/'
VAL_DIR = './validation/'

# Toggles
PLOT_TEMP = True
PLOT_FADEN_SCALARS = True
PLOT_FADEN_ITF = True

# Simulation Parameters
T_DELAY = 61.23  # s (Start delay due to initial LF)
T_DELAY_ITF = 61.2  # s (Slightly rounded delay used for interface extraction)
DT = 0.1  # s

AREA = 0.04 * 0.08  # m^2
AREA_2D = 0.04 * 1.0  # m^2
FULL_VOL = 0.04 * 0.04  # m^3
MASS_INI = 1.3886624  # kg


# =============================================================================
# 2. HELPER FUNCTIONS
# =============================================================================
def plot_with_faden(time_sim, y_sim, label_sim, time_fl, y_fl, label_fl,
                    faden_files=None, faden_labels=None, line_styles=None,
                    ylabel='', xlabel='Time [s]', scale_time=False,
                    exp_range_files=None, y_min=None, y_max=None):
    """Plots simulation, Fluent, and Faden data with discrete error bars."""
    plt.plot(time_sim, y_sim, '-k', label=label_sim)
    plt.plot(time_fl, y_fl, '--r', label=label_fl)

    # Plot specific Faden numerical/experimental points
    if PLOT_FADEN_SCALARS and faden_files:
        for file, style, label in zip(faden_files, line_styles, faden_labels):
            t_f, y_f = np.loadtxt(f'{VAL_DIR}{file}', skiprows=1, delimiter=',', unpack=True)
            if scale_time: t_f *= 60
            plt.plot(t_f, y_f, style, label=label)

    # Plot Faden experimental ranges as discrete error bars
    if PLOT_FADEN_SCALARS and exp_range_files:
        for file, label in exp_range_files:
            time_exp, LF_exp = np.loadtxt(f'{VAL_DIR}{file}', delimiter=',', skiprows=1, unpack=True)
            unique_times = np.unique(time_exp)

            LF_min = np.array([np.min(LF_exp[time_exp == t]) for t in unique_times])
            LF_max = np.array([np.max(LF_exp[time_exp == t]) for t in unique_times])

            LF_mean = (LF_max + LF_min) / 2
            LF_err = np.vstack((LF_mean - LF_min, LF_max - LF_mean))
            LF_err[LF_err == 0] = 0.002  # Small offset for visibility

            if scale_time: unique_times *= 60
            plt.errorbar(unique_times, LF_mean, yerr=LF_err, fmt='g.', capsize=3, label=label)

    plt.xlabel(xlabel, fontsize=16)
    plt.ylabel(ylabel, fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    if y_min is not None or y_max is not None:
        plt.ylim(y_min, y_max)
    plt.legend(fontsize=14)
    plt.tight_layout()
    plt.show()
    plt.close()


def extract_fluent_itf(filepath):
    """Efficiently reads Fluent interface .xy files by skipping headers/footers."""
    # Skip first 4 header lines and last footer line
    data = np.genfromtxt(filepath, skip_header=4, skip_footer=1)
    return data[:, 1], data[:, 0]  # returns x, y


def process_faden_itf(filepath, window=11, poly=2):
    """Loads and smooths the Faden interface data."""
    x, y = np.loadtxt(filepath, skiprows=1, delimiter=',', unpack=True)
    sorted_idx = np.argsort(y)
    x_smooth = savgol_filter(x[sorted_idx], window, poly)
    return x_smooth / 1000, y[sorted_idx] / 1000


# =============================================================================
# 3. SCALARS & ENERGY BALANCE
# =============================================================================
# Load report files
data_solid = np.loadtxt(f'{SIM_DIR}CFD_1/report-file.out', skiprows=3)
data_liquid = np.loadtxt(f'{SIM_DIR}CFD_2/report-file.out', skiprows=3)
data_fluent = np.loadtxt(f'{VAL_DIR}report-file.out', skiprows=3)

# Align arrays to the shortest simulation time
i_end = min(len(data_solid), len(data_liquid))
time = data_solid[:i_end, 7] + T_DELAY
time_l = data_liquid[:i_end, 9] + T_DELAY

# CoCoNuT variables
q_cool = -AREA * data_solid[:i_end, 1]
temp_s = data_solid[:i_end, 2:7]
vol_liquid = data_liquid[:i_end, 1]
v_max = data_liquid[:i_end, 2]
q_heat = AREA * data_liquid[:i_end, 3]
temp_l = data_liquid[:i_end, 4:9]

LF = vol_liquid / FULL_VOL
temp_combined = np.where(temp_l >= 311.13 - 3, np.minimum(temp_s, temp_l), np.maximum(temp_s, temp_l))

# Fluent variables
time_fl = data_fluent[:, 11]
LF_1_fl = 1 - (data_fluent[:, 1] / MASS_INI)
v_max_fl = data_fluent[:, 3]
q_heat_fl = AREA * data_fluent[:, 4]
q_cool_fl = -AREA * data_fluent[:, 5]
temp_fl = data_fluent[:, 6:11]

# Energy Balance
net_heat_flux = AREA_2D * (data_fluent[:, 4] + data_fluent[:, 5])
heat_flux_int = np.trapz(net_heat_flux, x=time_fl)

mass_end, h_hot = 1.3157895, 265038.54
dH = 149636.7 - (-18108.53)
H_out = (MASS_INI - mass_end) * h_hot
heat_net_in = heat_flux_int - H_out

print("\n--- ENERGY BALANCE ---")
print(f"Total heat transfer at walls: {heat_flux_int:.3f} J")
print(f"Total energy lost at outlet:  {H_out:.3f} J")
print(f"Net added energy:             {heat_net_in:.3f} J")
print(f"Enthalpy increase:            {dH:.3f} J")
print(f"Energy deficit:               {(abs(heat_net_in - dH) / dH) * 100:.3f} %\n")

# --- SCALAR PLOTTING ---
ts_0, ts_0_fl = math.ceil(1.0 / DT), 118

plot_with_faden(time[ts_0:], LF[ts_0:], 'Partitioned', time_fl[ts_0_fl:], LF_1_fl[ts_0_fl:], 'Fixed grid',
                faden_files=['Faden_LF_sim.csv'], faden_labels=['Faden - num.'], line_styles=['b.'],
                ylabel='Liquid fraction [-]', scale_time=True,
                exp_range_files=[('Faden_LF_exp_itv.csv', 'Faden - exp.')])

plot_with_faden(time[ts_0:], q_heat[ts_0:], 'Partitioned', time_fl[ts_0_fl:], q_heat_fl[ts_0_fl:], 'Fixed grid',
                faden_files=['Faden_HF_heated_exp.csv', 'Faden_HF_heated_sim.csv'],
                faden_labels=['Faden - exp.', 'Faden - num.'], line_styles=['g--', 'b--'],
                ylabel='Heat transfer rate [W]', y_min=2, y_max=4)

plot_with_faden(time[ts_0:], q_cool[ts_0:], 'Partitioned', time_fl[ts_0_fl:], q_cool_fl[ts_0_fl:], 'Fixed grid',
                faden_files=['Faden_HF_cooled_exp.csv', 'Faden_HF_cooled_sim.csv'],
                faden_labels=['Faden - exp.', 'Faden - num.'], line_styles=['g--', 'b--'],
                ylabel='Heat transfer rate [W]')

# Velocity Plot
plt.plot(time[ts_0:], v_max[ts_0:] * 1000, '-k', label='Partitioned')
plt.plot(time_fl[ts_0_fl:], v_max_fl[ts_0_fl:] * 1000, '--r', label='Fixed grid')
plt.ylabel('Max. velocity mag. [mm/s]', fontsize=16)
plt.xlabel('Time [s]', fontsize=16)
plt.xticks(fontsize=14);
plt.yticks(fontsize=14)
plt.legend(fontsize=16)
plt.tight_layout()
plt.show()
plt.close()

# Temperature Plot
if PLOT_TEMP:
    labels = ['u1', 'u2', 'u3', 'l1', 'l2']
    for i, lbl in enumerate(labels):
        plt.plot(time[ts_0:], temp_combined[ts_0:, i], '-k', label='Coconut')
        plt.plot(time_fl[ts_0_fl:], temp_fl[ts_0_fl:, i], '--r', label='Fluent')
        plt.ylabel(f'Temperature ({lbl}) [K]')
        plt.xlabel('Time [s]')
        plt.legend()
        plt.tight_layout()
        plt.show()
        plt.close()

# =============================================================================
# 4. INTERFACE TRACKING
# =============================================================================
print("--- INTERFACE PROCESSING ---")
pp = PostProcess(f'{SIM_DIR}case_results.pickle')
sx = pp.add_subset(interface='interface_x', model_part='boundary_in_nodes')
sy = pp.add_subset(interface='interface_y', model_part='boundary_out_faces')

x_coco = sx.get_values('coordinates', 'x')
y_coco = sx.get_values('coordinates', 'y')

itf_files_fl = ['itf-pos-257-96s.xy', 'itf-pos-738-97s.xy', 'itf-pos-1200-00s.xy', 'itf-pos-1800-00s.xy']

for itf in itf_files_fl:
    parts = itf.split('-')
    legend_time = f"{parts[2]}.{parts[3].replace('s.xy', '')}"
    t_target = round(float(parts[2]) + float(parts[3].replace('s.xy', '')) / 100, 1)

    # 1. CoCoNuT Partitioned Interface
    try:
        idx = int((t_target - T_DELAY_ITF) / DT)
        plt.plot(x_coco[idx, :].flatten(), y_coco[idx, :].flatten(), 'k-', label=f'Partitioned - {legend_time} s')
    except IndexError:
        print(f"Warning: CoCoNuT has not reached simulation time {t_target} s yet.")

    # 2. Fluent Enthalpy-Porosity Interface
    x_fl, y_fl = extract_fluent_itf(f'{VAL_DIR}{itf}')
    plt.plot(x_fl, y_fl, 'r--', label=f'Enthalpy-Porosity - {legend_time} s')

# 3. Faden Paper Interface
if PLOT_FADEN_ITF:
    x_fad, y_fad = process_faden_itf(f'{VAL_DIR}Faden-num-itf-1800s.csv')
    plt.plot(x_fad, y_fad, 'g--', label='Faden - 1800 s')

plt.ylabel('y-coordinate [m]', fontsize=16)
plt.xlabel('x-coordinate [m]', fontsize=16)
plt.xlim((0, 0.014))
plt.xticks(fontsize=14);
plt.yticks(fontsize=14)
# plt.legend(fontsize=16) # Left commented out as in original script
plt.tight_layout()
plt.show()
plt.close()