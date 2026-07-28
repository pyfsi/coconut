import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import medfilt
import os

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
PLOT_ANALYTICAL_STANDALONE = False  # Toggle for isolated analytical plots

# Hirata Specific Parameters
SIMULATIONS = [
    ("fluent_python", "Partitioned simulation", "-")
]
V_INITIAL = 0.052 * 0.052
VF_INITIAL = 0
TIME_OFFSET = 62.3

WALL_AREAS = {
    "htr-side": 0.104,  # m²
    "htr-bottom": 0.052,  # m²
    "htr-top": 0.052,  # m²
    "ht-int": 1.0  # Reference area
}

VALIDATION_FILES = [
    ("melt_frac_Hirata_sim.csv", "Kozak & Ziskind", "k--", "melt_fraction"),
    ("melt_frac_Hirata_exp.csv", "Hirata et al. (exp.)", "ro", "melt_fraction"),
    ("melt_frac_Hirata_ana.csv", "Analytical solution", "g-.", "melt_fraction"),
    ("sinking_velocity_hirata.csv", "Analytical solution", "k--", "v_y")
]

VAR_LABELS = {
    "vmax": "Maximum flow velocity [m/s]",
    "htr-bottom": "Heat transfer rate bottom [W]",
    "ht-int": "Heat transfer rate interface [W]",
    "htr-side": "Heat transfer rate side [W]",
    "htr-top": "Heat transfer rate top [W]",
    "v_y": "Solid y-velocity [m/s]",
    "melt_fraction": "Melt Fraction (Mass) [-]"
}

STYLE_MAP = {
    'Square case low tol': {'color': '#0033CC', 'linestyle': '-'},
    'Partitioned simulation': {'color': '#1F77B4', 'linestyle': '-'},
    'Kozak & Ziskind': {'color': '#000000', 'linestyle': '--'},
    'Hirata et al. (exp.)': {'color': '#CC0000', 'marker': 'o', 'linestyle': 'None'},
    'Analytical solution': {'color': '#008000', 'linestyle': '-.', 'linewidth': 1.5}
}

# Databases
rb_database = {}
rep_database = {}
validation_data = {}


# ===========================================================================
# 2. ANALYTICAL SOLUTION GENERATION
# ===========================================================================
def generate_analytical_data():
    """Generates analytical solution for sinking velocity and melt fraction."""
    H = 0.052
    g = 9.81
    rho_s = 856.0  # Solid density
    rho_l = 775.13  # Liquid density
    mu = 0.00401
    B = 0.052
    k = 0.15
    T_w = 329.85
    T_m = 309.85
    L = 248000.0

    t = np.arange(1, 601)  # Time from 1 to 600s

    term1 = (g * (rho_s - rho_l)) / (mu * B ** 2)
    term2 = (k * (T_w - T_m)) / (rho_s * L)
    term3 = rho_l / rho_s

    A = (term1 * (term2 ** 3) * term3) ** 0.25
    inner_term = H ** 0.75 - 0.75 * A * t

    # Calculate u_s (clip inner term to 0 to avoid complex numbers if it melts completely)
    u_s = np.zeros_like(t, dtype=float)
    valid = inner_term > 0
    u_s[valid] = A * (inner_term[valid]) ** (1 / 3)

    # Calculate s (position) and melt fraction (s/H)
    s = np.zeros_like(t, dtype=float)
    s[valid] = H - (inner_term[valid]) ** (4 / 3)
    s[~valid] = H  # If inner_term <= 0, it means the solid has completely melted
    melt_fraction = s / H

    if PLOT_ANALYTICAL_STANDALONE:
        # Plotting Sinking Velocity
        plt.figure(figsize=(10, 6))
        plt.plot(t, u_s * 1000, color='b', linewidth=2)
        plt.xlabel("Time (s)", fontsize=12)
        plt.ylabel("Sinking Velocity $u_s$ (mm/s)", fontsize=12)
        plt.title("Sinking Velocity vs. Time", fontsize=14)
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.tight_layout()
        plt.show()

        # Plotting Melt Fraction
        plt.figure(figsize=(10, 6))
        plt.plot(t, melt_fraction, color='r', linewidth=2)
        plt.xlabel("Time (s)", fontsize=12)
        plt.ylabel("Melt Fraction ($s/H$)", fontsize=12)
        plt.title("Melt Fraction vs. Time", fontsize=14)
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.tight_layout()
        plt.show()

    # CSV Writing
    os.makedirs(VALIDATION_DIR, exist_ok=True)

    with open(os.path.join(VALIDATION_DIR, "sinking_velocity_hirata.csv"), "w", newline="") as f:
        for time_val, us_val in zip(t, u_s):
            time_str = str(float(time_val)).replace('.', ',')
            us_str = str(-us_val).replace('.', ',')  # Negative sign retained for v_y mapping
            f.write(f"{time_str}; {us_str}\n")

    with open(os.path.join(VALIDATION_DIR, "melt_frac_Hirata_ana.csv"), "w", newline="") as f:
        for time_val, mf_val in zip(t, melt_fraction):
            time_str = str(float(time_val)).replace('.', ',')
            mf_str = str(mf_val).replace('.', ',')
            f.write(f"{time_str}; {mf_str}\n")

    print("  -> Analytical CSVs generated successfully.")


# ===========================================================================
# 3. HELPER FUNCTIONS & DATA CLEANING
# ===========================================================================
def add_to_db(db, var, label, t, vals, fmt):
    db.setdefault(var, []).append((label, t, vals, fmt))


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

    # Hirata specific column extraction
    headers = ["vmax", "htr-side", "htr-bottom", "htr-top", "ht-int"] if data.shape[1] == 7 else [f"Var_{k}" for k in
                                                                                                  range(data.shape[
                                                                                                            1] - 2)]

    for i, var in enumerate(headers):
        add_to_db(rep_database, var, label, time_rep[2:], data[2:, i + 1], fmt)


# ===========================================================================
# 5. PLOTTING & EXECUTION LOGIC
# ===========================================================================
def finalize_and_save_plot(fig, ax, var_name, filename, y_label):
    ax.set_xlabel("Time [s]", fontsize=18)
    ax.set_ylabel(y_label, fontsize=18)
    ax.tick_params(axis='both', labelsize=14)
    ax.grid(True, linestyle='--', alpha=0.7)

    lines, labels = ax.get_legend_handles_labels()
    fig.tight_layout()
    ax.legend(lines, labels, loc='best', fancybox=True, fontsize=14)

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
        fig, ax = plt.subplots(figsize=(8.0, 6.0))

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


def plot_htr_comparison(db, prefix=""):
    htr_vars = [v for v in db if v.startswith("htr-")]
    for sim in {l for var in htr_vars for l, _, _, _ in db[var]}:
        fig, ax = plt.subplots(figsize=(8.0, 6.0))
        for var in htr_vars:
            for label, t, vals, _ in db[var]:
                if label == sim:
                    is_flux = PLOT_HEAT_FLUX and var in WALL_AREAS
                    plot_vals = vals / WALL_AREAS[var] if is_flux else vals
                    line_label = VAR_LABELS.get(var, var).replace("Heat transfer rate ", "").replace(" [W]",
                                                                                                     "").capitalize()
                    if is_flux: line_label += " Flux"
                    ax.plot(t, plot_vals, label=line_label, linewidth=1.5)

        fname = f"{prefix}HTR_{sim.replace(' ', '_').replace('/', '_').replace('°', 'deg')}.png"
        y_axis_label = "Heat Flux [W/m²]" if PLOT_HEAT_FLUX else "Heat Transfer Rate [W]"
        finalize_and_save_plot(fig, ax, "htr", fname, y_axis_label)


# *** MAIN EXECUTION ***
if __name__ == "__main__":
    print("=== Running Hirata Validation Post-Processing ===")
    generate_analytical_data()
    load_validation_data()

    for folder, label, fmt in SIMULATIONS:
        print(f"*** Processing: {label} ({folder}) ***")
        process_rigid_body(folder, label, fmt)
        process_report_file(folder, label, fmt)

    print("\nGenerating Plots & Calculations...")
    print_time_averaged_htr(rep_database)
    plot_standard(rb_database, prefix="Solid")

    # Hirata case specific: Group the heat transfer rates together
    plot_standard(rep_database, prefix="", exclude=("htr-",))
    plot_htr_comparison(rep_database, prefix="")

    print("\nProcessing Complete.")
    if SHOW_PLOTS: plt.show()