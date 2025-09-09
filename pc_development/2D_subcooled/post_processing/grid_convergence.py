import numpy as np

# ======================
# Settings
# ======================
t_delay_001, t_delay_002 = 12.47, 61.23
start_lf_001 = False
t_delay = t_delay_001 if start_lf_001 else t_delay_002

# Directories
fine_case = '../Faden_split_3_rerun/run_6/'
coarse_case = '../Faden_split_3_coarse/run_7/'
report_file = 'report-file.out'
faden_dir = 'Faden_paper/'

# Geometry
area = 0.04 * 0.08
full_vol = 0.04 * 0.04

# Grid spacing (characteristic sizes)
h_fine = 1/800  # replace with actual fine mesh spacing
h_coarse = 1/400  # replace with actual coarse mesh spacing
r = h_coarse / h_fine  # refinement ratio


# ======================
# Helper functions
# ======================
def load_case(case_dir):
    solid = np.loadtxt(f'{case_dir}CFD_1/{report_file}', delimiter=' ', skiprows=3)
    liquid = np.loadtxt(f'{case_dir}CFD_2/{report_file}', delimiter=' ', skiprows=3)

    time = solid[:, 7] + t_delay
    time_l = liquid[:, 9] + t_delay
    n = min(len(time), len(time_l))

    q_cool = -area * solid[:n, 1]
    q_heat = area * liquid[:n, 3]
    vol_liquid = liquid[:n, 1]
    LF = vol_liquid / full_vol

    return time[:n], LF, q_heat, q_cool


def load_faden(filename, scale_time=True):
    data = np.loadtxt(f'{faden_dir}{filename}', delimiter=',', skiprows=1)
    t, y = data[:, 0], data[:, 1]
    if scale_time: t *= 60
    return t, y


def discrete_l2_error(ref_time, ref_values, time, values):
    interp = np.interp(ref_time, time, values)
    err = ref_values - interp
    return np.sqrt(np.sum(err ** 2) / np.sum(ref_values ** 2))


def roache_gci(fine_err, coarse_err, r, safety=1.25):
    """
    Compute Roache's Grid Convergence Index (GCI) using the generalized Richardson extrapolation.

    fine_err: error of fine mesh vs reference
    coarse_err: error of coarse mesh vs reference
    r: refinement ratio (coarse/fine)
    safety: safety factor (default 1.25)
    """
    # Avoid division by zero
    if fine_err == coarse_err:
        return 0.0, np.nan  # GCI, order undefined

    # Estimated order of convergence p
    p = np.log(abs(coarse_err / fine_err)) / np.log(r)

    # GCI on fine grid
    gci_fine = safety * fine_err / (r ** p - 1) * 100  # percentage
    return gci_fine, p


# ======================
# Load data
# ======================
time_fine, LF_fine, qh_fine, qc_fine = load_case(fine_case)
time_coarse, LF_coarse, qh_coarse, qc_coarse = load_case(coarse_case)

t_faden_LF, LF_faden = load_faden('Faden_LF_sim.csv')
t_faden_qh, qh_faden = load_faden('Faden_HF_heated_sim.csv')
t_faden_qc, qc_faden = load_faden('Faden_HF_cooled_sim.csv')

# ======================
# Discrete L2 Errors vs Faden
# ======================
lf_err_fine = discrete_l2_error(t_faden_LF, LF_faden, time_fine, LF_fine)
lf_err_coarse = discrete_l2_error(t_faden_LF, LF_faden, time_coarse, LF_coarse)

qh_err_fine = discrete_l2_error(t_faden_qh, qh_faden, time_fine, qh_fine)
qh_err_coarse = discrete_l2_error(t_faden_qh, qh_faden, time_coarse, qh_coarse)

qc_err_fine = discrete_l2_error(t_faden_qc, qc_faden, time_fine, qc_fine)
qc_err_coarse = discrete_l2_error(t_faden_qc, qc_faden, time_coarse, qc_fine)  # note: qc_fine typo fixed

# ======================
# Roache GCI
# ======================
lf_gci, lf_order = roache_gci(lf_err_fine, lf_err_coarse, r)
qh_gci, qh_order = roache_gci(qh_err_fine, qh_err_coarse, r)
qc_gci, qc_order = roache_gci(qc_err_fine, qc_err_coarse, r)

# ======================
# Report
# ======================
print("=== Discrete L2 Errors vs Faden ===")
print(f"LF: Fine {lf_err_fine:.6f}, Coarse {lf_err_coarse:.6f}")
print(f"Heat Flux Heated: Fine {qh_err_fine:.6f}, Coarse {qh_err_coarse:.6f}")
print(f"Heat Flux Cooled: Fine {qc_err_fine:.6f}, Coarse {qc_err_coarse:.6f}")

print("\n=== Roache GCI [%] and Estimated Order of Convergence ===")
print(f"LF: GCI {lf_gci:.2f}%, p = {lf_order:.2f}")
print(f"Heat Flux Heated: GCI {qh_gci:.2f}%, p = {qh_order:.2f}")
print(f"Heat Flux Cooled: GCI {qc_gci:.2f}%, p = {qc_order:.2f}")
