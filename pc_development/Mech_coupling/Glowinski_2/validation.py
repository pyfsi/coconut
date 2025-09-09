import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ----------------------
# File paths
# ----------------------
GLOW_YCO_PATH = "./post_processing/glowinski_data/y_co.csv"
GLOW_YVEL_PATH = "./post_processing/glowinski_data/y_vel.csv"
FLUENT_PATH = "./post_processing/new-report-file.out"
SIM_PATH = "./CFD_2/report-file.out"
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

def load_simulation_data(file_path):
    data = np.loadtxt(file_path, delimiter=' ', skiprows=3)

    # Extract columns (0-based indexing)
    y = data[:, 4]  # pos_y
    vel_y = data[:, 6]  # vel_y
    time = data[:, 7]  # flow_time

    return time, y, vel_y

# ----------------------
# Prepare output folder
# ----------------------
os.makedirs(FIG_DIR, exist_ok=True)

# ----------------------
# Load datasets
# ----------------------
glow_time, glow_y, glow_time_vel, glow_vel = load_glowinski_data(GLOW_YCO_PATH, GLOW_YVEL_PATH)
fluent_time, fluent_y, fluent_vel = load_simulation_data(FLUENT_PATH)
sim_time, sim_y, sim_vel = load_simulation_data(SIM_PATH)

# ----------------------
# Plot 1: Y vs Time
# ----------------------
plt.figure()
plt.plot(glow_time, glow_y/10, label='Glowinski', linestyle='--')
plt.plot(fluent_time, fluent_y*100, label='Fluent', linestyle='-')
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
plt.plot(sim_time, sim_vel*100, label='Partitioned', linestyle='-')
plt.xlabel('Time [s]')
plt.ylabel('Y-Velocity [cm/s]')
plt.legend()
plt.grid(True)

plt.savefig(os.path.join(FIG_DIR, "vel_vs_time.png"))
plt.show()