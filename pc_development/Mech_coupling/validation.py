import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter

# ----------------------
# Smoother
# ----------------------
smoothing_on = False
def smooth_signal(signal, window=51, polyorder=3):
    """
    Smooth a 1D signal using Savitzky-Golay filter.
    - window: must be odd, controls smoothing span
    - polyorder: polynomial order used for fitting
    """
    # Ensure window length is valid
    if window >= len(signal):
        window = len(signal) - 1 if len(signal) % 2 == 0 else len(signal)
    if window < 3:  # avoid too small windows
        return signal
    if window % 2 == 0:
        window += 1

    return savgol_filter(signal, window_length=window, polyorder=polyorder)

# ----------------------
# File paths
# ----------------------
# Validation data
GLOW_YCO_PATH = "./post_processing/glowinski_data/y_co.csv"
GLOW_YVEL_PATH = "./post_processing/glowinski_data/y_vel.csv"

# Fluent Simulation
FLUENT_PATH = "/cfdfile1/data/fm/victor/Documents/PC_wp2_dev/Terminal_1/report-file.out" # "./post_processing/report-file-fluent.out"
FLUENT_RB_PATH = "/cfdfile1/data/fm/victor/Documents/PC_wp2_dev/Terminal_1/disk_motion_disk.6dof"

# Partitioned simulation
SIM_PATH = "./Terminal_1/CFD_2/rigid-body-report-file.out"
SIM_FLUENT_PATH = "./Terminal_1/CFD_2/report-file.out"

FIG_DIR = "./post_processing/figures"

# ----------------------
# Load data
# ----------------------
def load_glowinski_data(yco_path, yvel_path):
    # Read raw data
    yco = pd.read_csv(yco_path, sep=';', header=None, names=['time', 'y'])
    yvel = pd.read_csv(yvel_path, sep=';', header=None, names=['time', 'vel_y'])

    # Convert comma decimals to floats
    yco['time'] = yco['time'].str.replace(',', '.').astype(float)
    yco['y'] = yco['y'].str.replace(',', '.').astype(float)
    yvel['time'] = yvel['time'].str.replace(',', '.').astype(float)
    yvel['vel_y'] = yvel['vel_y'].str.replace(',', '.').astype(float)

    # Return as numpy arrays
    return (
        yco['time'].to_numpy(),
        yco['y'].to_numpy(),
        yvel['time'].to_numpy(),
        yvel['vel_y'].to_numpy()
    )

def load_fluent_data(file_path):
    """Parser for Fluent-style report-file.out"""
    data = np.loadtxt(file_path, delimiter=' ', skiprows=3)

    # Extract columns (adjust according to Fluent output structure)
    y = data[:, 2]       # pos_y
    vel_y = data[:, 4]   # vel_y
    force_y = data[:, 5]  # force integral in y-dir.
    time = data[:, 6]    # flow_time

    return time, y, vel_y, force_y

def load_partitioned_data(file_path):
    """Parser for rigid-body-report-file.out created by update_report_file"""
    # Skip comment lines starting with "#"
    data = np.loadtxt(file_path, comments="#")

    time = data[:, 0]       # time
    y = data[:, 2]          # CG_Y
    vel_y = data[:, 4]      # V_Y
    force_y = data[:, 7]    # F_Y

    return time, y, vel_y, force_y

def load_fluent_rb_data(file_path):
    """Parser for *.6dof file created by Fluent"""
    # Skip comment lines starting with "#"
    data = np.loadtxt(file_path, comments="#")

    time = data[:, 0]    # time
    y = data[:, 2]       # CG_Y

    return time, y

# ----------------------
# Prepare output folder
# ----------------------
os.makedirs(FIG_DIR, exist_ok=True)

# ----------------------
# Load datasets
# ----------------------
glow_time, glow_y, glow_time_vel, glow_vel = load_glowinski_data(GLOW_YCO_PATH, GLOW_YVEL_PATH)
glow_y -= (glow_y[0] - 40) # Recalibrate the curve

fluent_time, fluent_y, fluent_vel, fluent_force = load_fluent_data(FLUENT_PATH) # Fluent
fluent_rb_time, fluent_rb_y = load_fluent_rb_data(FLUENT_RB_PATH)

sim_force_time, _, _, sim_force = load_fluent_data(SIM_FLUENT_PATH) # Partitioned
sim_time, sim_y, sim_vel, sim_force_coco = load_partitioned_data(SIM_PATH)
sim_y[0] = 0.04 # m, set initial condition, otherwise 0

if smoothing_on:
    sim_vel_smooth = smooth_signal(sim_vel, window=21, polyorder=3) # Smooth Partitioned velocity
else:
    sim_vel_smooth = sim_vel

# ----------------------
# Plot 1: Y vs Time
# ----------------------
plt.figure()
plt.plot(glow_time, glow_y/10, label='Glowinski', linestyle='--')
plt.plot(fluent_rb_time, fluent_rb_y*100, label='Fluent', linestyle='-')
plt.plot(sim_time, sim_y*100, label='Partitioned', linestyle='-')
plt.xlabel('Time [s]')
plt.ylabel('Y-coordinate [cm]')
plt.legend()
plt.grid(True)

plt.savefig(os.path.join(FIG_DIR, "y_vs_time.png"))
plt.show()

# ----------------------
# Plot 2: Y-Velocity vs Time
# ----------------------
plt.figure()
plt.plot(glow_time_vel, glow_vel, label='Glowinski', linestyle='--')
plt.plot(fluent_time, fluent_vel*100, label='Fluent', linestyle='-')
plt.plot(sim_time, sim_vel_smooth*100, label='Partitioned', linestyle='-')
plt.xlabel('Time [s]')
plt.ylabel('Y-Velocity [cm/s]')
plt.legend()
plt.grid(True)

plt.savefig(os.path.join(FIG_DIR, "vel_vs_time.png"))
plt.show()

# ----------------------
# Plot 3: Y-Force vs Time
# ----------------------
plt.figure()
plt.plot(fluent_time, fluent_force, label='Fluent', linestyle='-')
plt.plot(sim_time, sim_force_coco, label='Partitioned', linestyle='-')
plt.xlabel('Time [s]')
plt.ylabel('Y-Force integral [N]')
plt.legend()
plt.grid(True)

plt.savefig(os.path.join(FIG_DIR, "force_vs_time.png"))
plt.show()