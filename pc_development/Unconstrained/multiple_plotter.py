import numpy as np
import matplotlib.pyplot as plt
import os

# ===========================================================================
# CONFIGURATION
# ===========================================================================

# List of simulations to compare: (Folder Path, Legend Label)
simulations = [
    # ("sharp_volume/k9c1_2", "K9 C1.2"),
    # ("sharp_volume", "K9 C10 new disp"),
    # ("sharp_volume/k9c10_prev_disp", "K9 C10 old disp"),
    ("hpc/new_disp/new_k9_c1", "K9 C1"),
    ("hpc/new_disp/new_k9_c5", "K9 C5"),
    ("hpc/new_disp/new_k9_c10", "K9 C10"),
    ("hpc/new_disp/new_k9_c25", "K9 C25"),
    ("hpc/new_disp/new_k9_c50", "K9 C50"),
    ("hpc/new_disp/new_k9_c100", "K9 C100"),
    # ("Final_remesh", "Final Remesh"),
    # ("Final_flat_remesh", "Flat Remesh"),
    # ("Final_ext_remesh/k9c1_2", "Ext Remesh"),
]

# pre = "post_processing"  # Subfolder prefix
pre = "CFD_2"

# ===========================================================================
# DATA STORAGE
# ===========================================================================
# These dictionaries will hold data for all sims:
# Key = Variable Name, Value = List of tuples [(label, time_array, value_array), ...]
rb_database = {}
rep_database = {}


# ===========================================================================
# HELPER FUNCTIONS
# ===========================================================================

def is_near_zero(x, tol=1e-12):
    return abs(x) < tol


def add_to_database(database, var_name, label, time, values):
    if var_name not in database:
        database[var_name] = []
    database[var_name].append((label, time, values))


# ===========================================================================
# MAIN PROCESSING LOOP
# ===========================================================================

for folder, label in simulations:
    print(f"--- Processing: {label} ({folder}) ---")

    # -----------------------------------------------------------------------
    # 1. Process CoCoNuT Rigid Body Motion History
    # -----------------------------------------------------------------------
    filename_rb = os.path.join(folder, pre, "rigid-body-report-file.out")

    if os.path.exists(filename_rb):
        try:
            with open(filename_rb, "r") as f:
                lines = f.readlines()

            # 1.1 Extract Header
            header_line = None
            for line in lines:
                if line.strip().startswith("#") and "time" in line:
                    header_line = line
                    break

            labels_rb = []
            if header_line:
                labels_rb = header_line.replace("#", "").split()

            # 1.2 Find Data Start
            data_start = 0
            for i, line in enumerate(lines):
                stripped = line.strip()
                if stripped and not stripped.startswith("#") and stripped[0].isdigit():
                    data_start = i
                    break

            # 1.3 Load Data
            data_rb = np.loadtxt(filename_rb, skiprows=data_start)

            # 1.4 Remove "Fake Zero" Row
            remove_first_row = False
            if data_rb.ndim > 1:
                for col in range(1, data_rb.shape[1]):
                    first_val = data_rb[0, col]
                    second_val = data_rb[1, col]
                    if is_near_zero(first_val) and not is_near_zero(second_val):
                        remove_first_row = True
                        break

            if remove_first_row:
                data_rb = data_rb[1:, :]

            # 1.5 Store Data for Plotting
            time_rb = data_rb[:, 0]
            num_vars_rb = data_rb.shape[1]

            # Loop through columns (skipping time at index 0)
            for i in range(1, num_vars_rb):
                var_name = labels_rb[i] if len(labels_rb) > i else f"Var {i}"
                add_to_database(rb_database, var_name, label, time_rb, data_rb[:, i])

        except Exception as e:
            print(f"  Error processing rigid body file: {e}")
    else:
        print(f"  Warning: File not found {filename_rb}")

    # -----------------------------------------------------------------------
    # 2. Process Fluent Report File
    # -----------------------------------------------------------------------
    filename_rep = os.path.join(folder, pre, "report-file.out")

    if os.path.exists(filename_rep):
        try:
            # 2.1 Load Data
            data_rep = np.loadtxt(filename_rep, delimiter=' ', skiprows=3)

            # 2.2 Detect and Interpolate Zeros
            raw_time = data_rep[:, -1]
            middle_columns = data_rep[:, 1:-1]
            zeros_count = np.sum(middle_columns == 0.0, axis=1)
            bad_rows_mask = (zeros_count >= 2) & (raw_time > 0)

            if np.sum(bad_rows_mask) > 0:
                print(f"  -> Interpolating artifacts in {filename_rep}")
                valid_mask = ~bad_rows_mask
                valid_time = raw_time[valid_mask]

                for col in range(data_rep.shape[1]):
                    if col == data_rep.shape[1] - 1: continue  # Skip time
                    valid_data = data_rep[valid_mask, col]
                    data_rep[bad_rows_mask, col] = np.interp(
                        raw_time[bad_rows_mask], valid_time, valid_data
                    )

            time_rep = data_rep[:, -1]

            # 2.3 Determine Header
            cols = np.shape(data_rep)[1]
            if cols == 7:
                header = ["vmax", "htr-vert", "htr-bottom", "htr-arc", "htr-int"]
            elif cols == 8:
                header = ["vmax", "htr-vert", "htr-bottom", "htr-arc", "htr-int", "vol-air"]
            elif cols == 9:
                header = ["vmax", "htr-vert", "htr-bottom", "htr-arc", "htr-int", "vol-air", "mass-source"]
            else:
                header = [f"Var_{k}" for k in range(cols - 2)]

            # 2.4 Store Data for Plotting
            # Note: Original script plotted from index 2 onwards (data_rep[2:, ...])
            start_idx = 2
            for i in range(len(header)):
                var_name = header[i]
                # Column index in data_rep is i+1 (because col 0 is iterations)
                val_data = data_rep[start_idx:, i + 1]
                t_data = time_rep[start_idx:]
                add_to_database(rep_database, var_name, label, t_data, val_data)

        except Exception as e:
            print(f"  Error processing report file: {e}")
    else:
        print(f"  Warning: File not found {filename_rep}")

# ===========================================================================
# PLOTTING
# ===========================================================================
print("\nGenerating Plots...")


# Function to plot a database dictionary
def plot_all_variables(database, prefix=""):
    for var_name, data_list in database.items():
        plt.figure()

        # Plot each simulation line
        for label, time, values in data_list:
            plt.plot(time, values, label=label)  # removed hardcoded linestyle/color

        plt.xlabel("Time [s]")
        plt.ylabel(var_name)
        plt.title(f"{prefix}: {var_name}")
        plt.grid(True)
        plt.legend()  # Show labels


# 1. Plot Rigid Body Variables
plot_all_variables(rb_database, prefix="RB")

# 2. Plot Report File Variables
plot_all_variables(rep_database, prefix="Report")

plt.show()