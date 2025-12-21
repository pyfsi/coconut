import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import trapezoid

folder = "Full_contact_f"

# ===========================================================================
# PART 1: Process CoCoNuT Rigid Body Motion History
# ===========================================================================

filename_rb = folder + "/CFD_2/rigid-body-report-file.out"
print(f"Processing: {filename_rb}")

try:
    with open(filename_rb, "r") as f:
        lines = f.readlines()

    # ---- 1.1 Extract Header ----------------------------------------------
    header_line = None
    for line in lines:
        if line.strip().startswith("#") and "time" in line:
            header_line = line
            break

    if header_line is None:
        print(f"Warning: Could not find header in {filename_rb}")
        labels_rb = []
    else:
        labels_rb = header_line.replace("#", "").split()

    # ---- 1.2 Find Data Start ---------------------------------------------
    data_start = 0
    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped and not stripped.startswith("#") and stripped[0].isdigit():
            data_start = i
            break

    # ---- 1.3 Load Data ---------------------------------------------------
    data_rb = np.loadtxt(filename_rb, skiprows=data_start)


    # ---- 1.4 Remove "Fake Zero" Row --------------------------------------
    def is_near_zero(x, tol=1e-12):
        return abs(x) < tol


    remove_first_row = False
    # Check each column except time (col 0)
    if data_rb.ndim > 1:
        for col in range(1, data_rb.shape[1]):
            first_val = data_rb[0, col]
            second_val = data_rb[1, col]
            if is_near_zero(first_val) and not is_near_zero(second_val):
                remove_first_row = True
                break

    if remove_first_row:
        print(" -> Removing artificial zero timestep.")
        data_rb = data_rb[1:, :]

    # ---- 1.5 Plot Rigid Body Data ----------------------------------------
    time_rb = data_rb[:, 0]
    num_vars_rb = data_rb.shape[1]

    for i in range(1, num_vars_rb):
        plt.figure()
        plt.plot(time_rb, data_rb[:, i], label='Rigid Body Data')
        plt.xlabel(labels_rb[0] if labels_rb else "Time")
        plt.ylabel(labels_rb[i] if len(labels_rb) > i else f"Var {i}")
        plt.title(f"RB: {labels_rb[i]} vs {labels_rb[0]}")
        plt.grid(True)

        if labels_rb[i] == "volume":
            volume = data_rb[:, i].copy()

except Exception as e:
    print(f"Error processing rigid body file: {e}")

# ===========================================================================
# PART 2: Process Fluent Report File
# ===========================================================================

filename_rep = folder + "/CFD_2/report-file.out"
print(f"\nProcessing: {filename_rep}")

try:
    # ---- 2.2 Load Data ---------------------------------------------------
    data_rep = np.loadtxt(filename_rep, delimiter=' ', skiprows=3)
    time_rep = data_rep[:, 6]
    header = ["vmax", "htr-vert", "htr-bottom", "htr-arc", "htr-int"]

    # ---- 2.4 Plot Report Data --------------------------------------------
    htr = None

    for i in range(len(header)):
        plt.figure()
        plt.plot(time_rep, data_rep[:, i+1], color='tab:orange', linestyle='--')
        plt.xlabel("Time [s]")
        plt.ylabel(header[i])
        plt.grid(True)

        if header[i] == "htr-int":
            htr = data_rep[:, i+1].copy()

except Exception as e:
    print(f"Error processing report file: {e}")

# ===========================================================================
# PART 3: Check energy balance of solid
# ===========================================================================
if htr is not None:
    rho = 775.13 # kg/m^3
    L = 248000.0 # J/kg

    end = min(np.size(volume), np.size(htr)) - 1
    latent = rho * L * (volume[end] - volume[1])
    htr_int = trapezoid(htr[1:end], time_rep[1:end])
    print(f"Stored latent heat = {latent} J")
    print(f"Total heat transfer = {htr_int} J")
    print(f'Ratio Stored/Added = {latent/htr_int}')

# ===========================================================================
# Final Display
# ===========================================================================
plt.show()