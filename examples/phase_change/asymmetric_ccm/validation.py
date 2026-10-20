import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import medfilt
import os
from matplotlib.lines import Line2D

# ===========================================================================
# 1. DIRECTORIES & CONFIGURATION
# ===========================================================================
# Folder paths
PRE_FOLDER = "post_processing"  # Subfolder inside each simulation folder containing .out files
OUTPUT_DIR = os.path.join("validation", "figures")
VALIDATION_DIR = "validation"  # Folder containing the experimental CSV files

# Toggles
SAVE_PLOTS = False
SHOW_PLOTS = True
SHOW_VALIDATION = True
PRINT_TIME_AVERAGED_HTR = True
PLOT_HEAT_FLUX = True  # True plots [W/m²], False plots [W]

# 0° Mode Specific Parameters
SIMULATIONS = [
    ("fluent_python", "10W BC, var. density, 0.2 mm", "-")
]
V_INITIAL = 0.0028033717
VF_INITIAL = 0.014961
TIME_OFFSET = 25.532

WALL_AREAS = {
    "htr-vert": 0.050825585,  # m²
    "htr-bottom": 0.051800002,  # m²
    "htr-arc": 0.098426465,  # m²
    "ht-int": 1.0  # Reference area
}

VALIDATION_FILES = [
    ("melt_frac_exp_0.csv", "Schockner et al.", "k--", "melt_fraction")
]

VAR_LABELS = {
    "vmax": "Maximum flow velocity [m/s]",
    "htr-vert": "Heat transfer rate left [W]",
    "htr-bottom": "Heat transfer rate bottom [W]",
    "htr-arc": "Heat transfer rate arc [W]",
    "ht-int": "Heat transfer rate interface [W]",
    "vol-air": "Air volume [m³]",
    "mass-source": "Mass source [kg/s]",
    "time": "Time [s]",
    "cg_x": "CG x-coordinate [m]",
    "cg_y": "CG y-coordinate Y [m]",
    "theta_z": "Rotation Angle [°]",
    "v_x": "Solid x-velocity [m/s]",
    "v_y": "Solid y-velocity [m/s]",
    "f_x": "Fluid force (x) [N]",
    "f_y": "Fluid force (y) [N]",
    "m_z": "Fluid moment [Nm]",
    "melt_fraction": "Melt Fraction (Mass) [-]",
    "h_min": "Minimum gap width [m]"
}

STYLE_MAP = {
    'T BC, cst. density, 0.2 mm': {'color': '#0033CC', 'linestyle': '-'},
    'T BC, cst. density, 0.25 mm': {'color': '#809FFF', 'linestyle': '--'},
    'T BC, var. density, 0.2 mm': {'color': '#008000', 'linestyle': '-'},
    'T BC, var. density, 0.25 mm': {'color': '#80CC80', 'linestyle': '--'},
    '10W BC, var. density, 0.2 mm': {'color': '#CC0000', 'linestyle': '-'},
    '10W BC, var. density, 0.25 mm': {'color': '#FF8080', 'linestyle': '--'},
    'Schockner et al.': {'color': '#000000', 'linestyle': '-.', 'linewidth': 1.5}
}

# Databases
rb_database = {}
rep_database = {}
validation_data = {}


# ===========================================================================
# 2. HELPER FUNCTIONS & DATA CLEANING
# ===========================================================================
def add_to_db(db, var, label, t, vals, fmt):
    db.setdefault(var, []).append((label, t, vals, fmt))


def get_fluent_header(cols):
    headers = {
        7: ["vmax", "htr-vert", "htr-bottom", "htr-arc", "ht-int"],
        8: ["vmax", "htr-vert", "htr-bottom", "htr-arc", "ht-int", "vol-air"],
        9: ["vmax", "htr-vert", "htr-bottom", "htr-arc", "ht-int", "vol-air", "mass-source"]
    }
    return headers.get(cols, [f"Var_{k}" for k in range(cols - 2)])


def clean_remeshing_artifacts(y, kernel_size=0, threshold_sigma=3.0, recovery_window=8):
    if kernel_size % 2 == 0: kernel_size += 1
    residual = np.abs(y - medfilt(y, kernel_size))
    mad = np.median(residual[residual > 1e-12]) if np.median(residual) < 1e-12 and len(
        residual[residual > 1e-12]) > 0 else max(np.median(residual), 1e-12)

    spike_mask = residual > (threshold_sigma * 1.4826 * mad)
    bad_mask = np.copy(spike_mask)
    for idx in np.where(spike_mask)[0]:
        bad_mask[max(0, idx - 1):min(len(y), idx + recovery_window + 1)] = True

    y_clean = np.copy(y)
    if np.any(bad_mask) and np.sum(~bad_mask) > 2:
        x = np.arange(len(y))
        y_clean[bad_mask] = np.interp(x[bad_mask], x[~bad_mask], y[~bad_mask])
    return y_clean


def load_validation_data():
    if not SHOW_VALIDATION: return
    for f_name, label, fmt, var in VALIDATION_FILES:
        path = os.path.join(VALIDATION_DIR, f_name)
        if not os.path.exists(path): continue

        try:
            data = np.loadtxt(path, delimiter=chr(59), converters={
                0: lambda s: float(s.decode().replace(',', '.')),
                1: lambda s: float(s.decode().replace(',', '.'))
            })
            if data.size > 0:
                validation_data.setdefault(var, []).append((data[:, 0], data[:, 1], label, fmt))
        except Exception as e:
            print(f"  Error loading {path}: {e}")


# ===========================================================================
# 3. TRAJECTORY SPECIFIC FUNCTIONS
# ===========================================================================
def rotate_point(p, angle_deg):
    """Rotates a point (x, y) around the origin (0,0) by a given angle in degrees."""
    angle_rad = np.radians(angle_deg)
    x, y = p
    new_x = x * np.cos(angle_rad) - y * np.sin(angle_rad)
    new_y = x * np.sin(angle_rad) + y * np.cos(angle_rad)
    return (new_x, new_y)


def load_trajectory_validation_data(filename):
    """
    Loads validation CSV data for trajectories.
    Handles comma decimals and semicolon delimiters. Converts from mm to m.
    """
    if not os.path.exists(filename):
        print(f"Validation file not found: {filename}")
        return None

    try:
        with open(filename, 'r') as f:
            lines = f.readlines()

        data = []
        for line in lines:
            line = line.strip()
            if not line:
                continue
            parts = line.replace(',', '.').split(';')
            if len(parts) >= 2:
                data.append([float(parts[0]), float(parts[1])])

        return np.array(data) * 1e-3
    except Exception as e:
        print(f"Error parsing validation file {filename}: {e}")
        return None


def calculate_trajectories_with_domain(filename, p1_start, p2_start,
                                       domain_rotation_deg=0.0,
                                       enable_validation=False,
                                       enable_links=False,
                                       validation_dir=''):
    """
    Calculates trajectories using robust differentiation and Heun's method.
    Plots the trajectories against the specified, rotatable domain boundaries.
    """
    try:
        data = np.loadtxt(filename, comments='#')
    except Exception as e:
        print(f"Error reading file: {e}")
        return

    # Fix Corrupted Data (Extrapolation)
    valid_data = data[1:, :]
    row_1 = valid_data[0, :]
    row_2 = valid_data[1, :]
    dt_initial = row_2[0] - row_1[0]

    extrapolated_row = 2 * row_1 - row_2
    extrapolated_row[0] = row_1[0] - dt_initial

    data = np.vstack([extrapolated_row, valid_data])
    print(f"Data fixed. Extrapolated start time: t={data[0, 0]:.4f}s")

    # Extract columns
    t = data[:, 0]
    cg_x = data[:, 1]
    cg_y = data[:, 2]
    v_x = data[:, 3]
    v_y = data[:, 4]
    theta_deg = data[:, 5]

    theta_rad = np.radians(theta_deg)
    omega = np.gradient(theta_rad, t)

    # Setup Initial Coordinates & Validation
    val_left = None
    val_right = None

    if enable_validation:
        rot_int = int(domain_rotation_deg)
        val_left_file = os.path.join(validation_dir, f'traj_left_{rot_int}.csv')
        val_right_file = os.path.join(validation_dir, f'traj_right_{rot_int}.csv')

        val_left = load_trajectory_validation_data(val_left_file)
        val_right = load_trajectory_validation_data(val_right_file)

    p1_rot = rotate_point(p1_start, domain_rotation_deg)
    p2_rot = rotate_point(p2_start, domain_rotation_deg)

    if enable_validation and val_left is not None and val_right is not None:
        p1_rot = val_left[0]
        p2_rot = val_right[0]
        print(f"Start points overridden by validation data ({rot_int}° files).")

    # Integration Function (Heun's Method)
    def integrate_path_heun(start_x, start_y):
        path_x = np.zeros(len(t))
        path_y = np.zeros(len(t))
        path_x[0] = start_x
        path_y[0] = start_y

        for i in range(len(t) - 1):
            dt = t[i + 1] - t[i]

            r_x_curr = path_x[i] - cg_x[i]
            r_y_curr = path_y[i] - cg_y[i]

            v_rot_x_curr = -omega[i] * r_y_curr
            v_rot_y_curr = omega[i] * r_x_curr

            vel_x_curr = v_x[i] + v_rot_x_curr
            vel_y_curr = v_y[i] + v_rot_y_curr

            pred_x = path_x[i] + vel_x_curr * dt
            pred_y = path_y[i] + vel_y_curr * dt

            r_x_next = pred_x - cg_x[i + 1]
            r_y_next = pred_y - cg_y[i + 1]

            v_rot_x_next = -omega[i + 1] * r_y_next
            v_rot_y_next = omega[i + 1] * r_x_next

            vel_x_next = v_x[i + 1] + v_rot_x_next
            vel_y_next = v_y[i + 1] + v_rot_y_next

            path_x[i + 1] = path_x[i] + (vel_x_curr + vel_x_next) / 2.0 * dt
            path_y[i + 1] = path_y[i] + (vel_y_curr + vel_y_next) / 2.0 * dt

        return path_x, path_y

    p1_x, p1_y = integrate_path_heun(p1_rot[0], p1_rot[1])
    p2_x, p2_y = integrate_path_heun(p2_rot[0], p2_rot[1])

    # Plotting
    plt.figure(figsize=(10, 10))

    R_domain = 51.8e-3
    start_angle = np.radians(-30 + domain_rotation_deg)
    end_angle = np.radians(90 + domain_rotation_deg)
    theta_domain = np.linspace(start_angle, end_angle, 200)

    dom_x = R_domain * np.cos(theta_domain)
    dom_y = R_domain * np.sin(theta_domain)

    domain_x = np.concatenate(([0], dom_x, [0]))
    domain_y = np.concatenate(([0], dom_y, [0]))

    plt.plot(domain_x, domain_y, 'k-', linewidth=2.5)
    plt.fill(domain_x, domain_y, 'k', alpha=0.05)

    plt.plot(p1_x, p1_y, 'r-', linewidth=2)
    plt.plot(p2_x, p2_y, 'b-', linewidth=2)

    if enable_validation:
        if val_left is not None:
            plt.plot(val_left[:, 0], val_left[:, 1], color='maroon', linestyle=':', linewidth=2.5, alpha=0.8)
        if val_right is not None:
            plt.plot(val_right[:, 0], val_right[:, 1], color='navy', linestyle=':', linewidth=2.5, alpha=0.8)

    plt.scatter([p1_x[0], p2_x[0]], [p1_y[0], p2_y[0]], c='g', s=60, zorder=5)
    plt.scatter([p1_x[-1], p2_x[-1]], [p1_y[-1], p2_y[-1]], c='k', marker='x', s=60, zorder=5)

    if enable_links:
        indices_to_plot = [0]
        max_time_multiple = int(t[-1] // 100.0)
        target_times = [i * 100.0 for i in range(max_time_multiple + 1)]

        for target_t in target_times:
            if target_t >= t[0] - 1e-5 and target_t <= t[-1] + 1e-5:
                idx = np.argmin(np.abs(t - target_t))
                if idx not in indices_to_plot:
                    indices_to_plot.append(idx)

        for idx in indices_to_plot:
            plt.plot([p1_x[idx], p2_x[idx]], [p1_y[idx], p2_y[idx]],
                     color='purple', linestyle='--', linewidth=1.5, alpha=0.7)
            plt.scatter([p1_x[idx], p2_x[idx]], [p1_y[idx], p2_y[idx]],
                        color='purple', s=15, zorder=4)

    legend_elements = [
        Line2D([0], [0], color='black', linestyle='-', lw=2, label='Partitioned simulation'),
        Line2D([0], [0], color='black', linestyle=':', lw=2.5, label='Experimental - Shockner et al.')
    ]
    plt.legend(handles=legend_elements, loc='best', fontsize='medium')

    plt.xlabel('X Position (m)')
    plt.ylabel('Y Position (m)')

    display_name = os.path.basename(filename)
    plt.title(f'Trajectory within Domain\n{display_name}')
    plt.axis('equal')
    plt.grid(True, linestyle=':', alpha=0.6)
    plt.tight_layout()

    if SAVE_PLOTS:
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        plt.savefig(os.path.join(OUTPUT_DIR, f"Trajectory_{display_name.replace('.out', '')}.png"), dpi=300)

    if SHOW_PLOTS:
        plt.show()
    else:
        plt.close()


# ===========================================================================
# 4. FILE PROCESSING FUNCTIONS
# ===========================================================================
def process_rigid_body(folder, label, fmt):
    path = os.path.join(folder, PRE_FOLDER, "rigid-body-report-file.out")
    if not os.path.exists(path): return

    with open(path, "r") as f:
        lines = f.readlines()

    headers = next((l.replace("#", "").split() for l in lines if l.strip().startswith("#") and "time" in l), [])
    data = np.loadtxt(path,
                      skiprows=next(i for i, l in enumerate(lines) if l.strip() and not l.strip().startswith("#")))

    if data.ndim > 1 and abs(data[0, 1]) < 1e-12 and abs(data[1, 1]) >= 1e-12:
        data = data[1:]

    time_rb = data[:, 0] + TIME_OFFSET
    for i in range(1, data.shape[1]):
        var = headers[i].lower() if len(headers) > i else f"var_{i}"
        clean_val = clean_remeshing_artifacts(data[:, i])

        if var == "volume":
            melt_frac = 1.0 - (clean_val / V_INITIAL)
            if VF_INITIAL: melt_frac += VF_INITIAL - melt_frac[0]
            add_to_db(rb_database, "melt_fraction", label, time_rb, melt_frac, fmt)
        else:
            add_to_db(rb_database, var, label, time_rb, clean_val, fmt)


def process_report_file(folder, label, fmt):
    path = os.path.join(folder, PRE_FOLDER, "report-file.out")
    if not os.path.exists(path): return

    data = np.loadtxt(path, delimiter=' ', skiprows=3)
    raw_time = data[:, -1]

    bad_rows = (np.sum(data[:, 1:-1] == 0.0, axis=1) >= 2) & (raw_time > 0)
    for idx in np.where(bad_rows)[0]: bad_rows[idx:min(len(bad_rows), idx + 5)] = True

    for col in range(data.shape[1] - 1):
        if np.sum(~bad_rows) > 2:
            data[bad_rows, col] = np.interp(raw_time[bad_rows], raw_time[~bad_rows], data[~bad_rows, col])
        data[:, col] = clean_remeshing_artifacts(data[:, col])

    time_rep = raw_time + TIME_OFFSET
    headers = get_fluent_header(data.shape[1])

    for i, var in enumerate(headers):
        add_to_db(rep_database, var, label, time_rep[2:], data[2:, i + 1], fmt)


# ===========================================================================
# 5. PLOTTING LOGIC
# ===========================================================================
def finalize_and_save_plot(fig, ax, var_name, filename, y_label):
    ax.set_xlabel("Time [s]", fontsize=18)
    ax.set_ylabel(y_label, fontsize=18)
    ax.tick_params(axis='both', labelsize=14)
    ax.grid(True, linestyle='--', alpha=0.7)

    lines, labels = ax.get_legend_handles_labels()

    margin_height = 0.24
    fig.tight_layout(rect=[0, margin_height, 1, 1])
    fig.legend(lines, labels, loc='upper center', bbox_to_anchor=(0.5, margin_height),
               fancybox=True, ncol=3 if "HTR" in filename else 2, fontsize=14)

    if SAVE_PLOTS:
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        plt.savefig(os.path.join(OUTPUT_DIR, filename), dpi=300)

    if not SHOW_PLOTS:
        plt.close(fig)


def print_time_averaged_htr(db):
    if not PRINT_TIME_AVERAGED_HTR: return
    title = "Time-Averaged Heat Flux" if PLOT_HEAT_FLUX else "Time-Averaged Heat Transfer Rates"
    print(f"\n=== {title} ===")

    for sim in sorted({l for var in db if var.startswith("htr-") for l, _, _, _ in db[var]}):
        print(f"\nSimulation: {sim}")
        for var in sorted(v for v in db if v.startswith("htr-")):
            for l, t, vals, _ in db[var]:
                if l == sim:
                    is_flux = PLOT_HEAT_FLUX and var in WALL_AREAS
                    calc_vals = vals / WALL_AREAS[var] if is_flux else vals
                    avg = np.trapz(calc_vals, t) / (t[-1] - t[0]) if len(t) > 1 and (t[-1] - t[0]) > 0 else np.mean(
                        calc_vals)

                    base_name = VAR_LABELS.get(var, var).replace('Heat transfer rate ', '').replace(' [W]',
                                                                                                    '').capitalize()
                    print_name = f"{base_name} Flux" if is_flux else base_name
                    unit = "W/m²" if is_flux else "W"
                    print(f"  * {print_name}: {avg:.4f} {unit}")


def plot_standard(db, prefix="", exclude=()):
    for var, data_list in db.items():
        if any(var.startswith(ex) for ex in exclude): continue
        fig, ax = plt.subplots(figsize=(7.8, 6.8))

        is_ht_var = var.startswith("htr-") or var == "ht-int"
        is_flux = PLOT_HEAT_FLUX and is_ht_var and var in WALL_AREAS

        for label, time, vals, fmt in data_list:
            plot_vals = vals / WALL_AREAS[var] if is_flux else vals
            style = {'linewidth': 1.5, **STYLE_MAP.get(label, {})}
            if style:
                ax.plot(time, plot_vals, label=label, **style)
            else:
                ax.plot(time, plot_vals, fmt, label=label, linewidth=1.5)

        for vt, vy, vlbl, vfmt in validation_data.get(var, []):
            plot_vy = vy / WALL_AREAS[var] if is_flux else vy
            style = {'markersize': 4, **STYLE_MAP.get(vlbl, {})}
            if style:
                ax.plot(vt, plot_vy, label=vlbl, **style)
            else:
                ax.plot(vt, plot_vy, vfmt, label=vlbl, markersize=4)

        fname = f"{prefix}_{var.replace('/', '_').replace(' ', '_')}.png".strip('_')
        y_label = VAR_LABELS.get(var, var)
        if is_flux: y_label = y_label.replace("transfer rate", "flux").replace("[W]", "[W/m²]")

        finalize_and_save_plot(fig, ax, var, fname, y_label)


# ===========================================================================
# 6. MAIN EXECUTION
# ===========================================================================
if __name__ == "__main__":
    print("=== Running 0° Mode Post-Processing ===")
    load_validation_data()

    for folder, label, fmt in SIMULATIONS:
        print(f"*** Processing: {label} ({folder}) ***")
        process_rigid_body(folder, label, fmt)
        process_report_file(folder, label, fmt)

    print("\nGenerating Plots & Calculations...")
    print_time_averaged_htr(rep_database)

    plot_standard(rb_database, prefix="Solid")
    plot_standard(rep_database, prefix="")

    print("\n=== Calculating Trajectories ===")
    # Setup parameters for trajectory mapping
    trajectory_filename = 'fluent_python/post_processing/rigid-body-report-file.out'
    DOMAIN_ROTATION = 0.0
    ENABLE_TRAJ_VALIDATION = True
    ENABLE_LINKS = False
    point_1 = (0.03138, 0.02917)
    point_2 = (0.03715, 0.01578)

    if os.path.exists(trajectory_filename):
        calculate_trajectories_with_domain(
            trajectory_filename,
            point_1,
            point_2,
            domain_rotation_deg=DOMAIN_ROTATION,
            enable_validation=ENABLE_TRAJ_VALIDATION,
            enable_links=ENABLE_LINKS,
            validation_dir=VALIDATION_DIR
        )
    else:
        print(f"Trajectory plotting skipped. File not found: {trajectory_filename}")

    print("\nProcessing Complete.")
    if SHOW_PLOTS: plt.show()