import matplotlib.pyplot as plt
import numpy as np
import io
from scipy.signal import savgol_filter
from coconut.examples.post_processing.post_processing import PostProcess

# ---------------- Settings ----------------
run = 'split_4'  # 'coarse', 'rerun', 'auto' or 'split_4

# Initial delay due to starting with LF < 0.02
t_001 = 12.5  # s
t_002 = 61.2  # s
t_run_1 = t_002  # initial offset

# Simulation times per run
if run == 'coarse':
    sim_times = [1800, 1620, 900, 1440, 720, 540, 180]
elif run == 'rerun':
    sim_times = [1800, 1080, 1080, 1260, 1440, 5400]
elif run == 'auto':
    sim_times = [1440, 720, 720, 720, 720]
elif run == 'split_4':
    sim_times = [1440, 1260, 900, 1440]

if run == 'coarse':
    common_path = '../Faden_split_3_coarse/'
elif run == 'rerun':
    common_path = '../Faden_split_3_rerun/'
elif run == 'auto':
    common_path = '../Faden_auto_remesh/'
elif run == 'split_4':
    common_path = '../Faden_split_4/'

# Time instances to plot
t_sims = [3600.0]  # seconds

# Partitioned simulation settings
dt = 0.1

# Fluent interface settings
plot_fluent = True
common_path_fl = './fluent_interfaces/'
# build: "itf-pos-{time}-00s.xy"
itf_files_fl = [f"itf-pos-{int(t)}-00s.xy" for t in t_sims]

# Faden paper settings
plot_Faden = True
common_path_Fa = './Faden_paper/'
# build: "Faden-num-itf-{time}s.csv", then "Faden-exp-itf-{time}s.csv"
itf_files_Fa = [f"Faden-{kind}-itf-{int(t)}s.csv"
                for t in t_sims
                for kind in ["num", "exp"]]

line_style = lambda file: 'r--' if 'num' in file else 'b--'


# ---------------- Plotting ----------------
lines = []

# Precompute cumulative offsets including initial delay
cumulative_times = np.insert(np.cumsum(sim_times) + t_run_1, 0, t_run_1)

# ---------------- Partitioned simulation ----------------
for i, t_sim in enumerate(t_sims):
    run_index = np.searchsorted(cumulative_times, t_sim, side='right') - 1
    t_local = t_sim - cumulative_times[run_index]
    case_file = f'{common_path}run_{run_index+1}/case_results.pickle'

    try:
        pp = PostProcess(case_file)
    except FileNotFoundError:
        print(f"File not found: {case_file}. Skipping this time instance.")
        continue

    sx = pp.add_subset(interface='interface_x', model_part='boundary_in_nodes')
    sy = pp.add_subset(interface='interface_y', model_part='boundary_out_faces')

    x = sx.get_values('coordinates', 'x')
    y = sx.get_values('coordinates', 'y')

    try:
        idx = int(t_local / dt)
        x_slice = x[idx, :].flatten()
        y_slice = y[idx, :].flatten()
        line, = plt.plot(x_slice, y_slice, 'k-')
        lines.append(line)
    except IndexError:
        print(f"Partitioned run {run_index+1} has not reached {t_local:.1f} s (local time).")

# ---------------- Fluent interfaces ----------------
if plot_fluent:
    for i, itf_file in enumerate(itf_files_fl):
        with open(common_path_fl + itf_file, 'r') as f:
            data_lines = f.readlines()[4:-1]
        data = np.loadtxt(io.StringIO(''.join(data_lines)), delimiter='\t')
        y_fl, x_fl = data[:,0], data[:,1]
        line, = plt.plot(x_fl, y_fl, 'g--')
        lines.append(line)

# ---------------- Faden interfaces ----------------
if plot_Faden:
    for i, itf_file in enumerate(itf_files_Fa):
        x, y = np.loadtxt(common_path_Fa + itf_file, skiprows=1, delimiter=',', unpack=True)
        sorted_idx = np.argsort(y)
        x, y = x[sorted_idx], y[sorted_idx]
        x_smooth = savgol_filter(x, 11, 2)
        line, = plt.plot(x_smooth/1000, y/1000, line_style(itf_file))
        lines.append(line)

# ---------------- Final plot ----------------
plt.xlabel('x-coordinate [m]', fontsize=16)
plt.ylabel('y-coordinate [m]', fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlim(0, 0.04)
plt.ylim(0, 0.04)

# Create legend handles manually to match colors/styles
legend_handles = []
k = 0

# Partitioned (take the first Partitioned line as representative)
legend_handles.append(lines[0])
k += len(t_sims)

# Fluent (first Fluent line)
if plot_fluent:
    legend_handles.append(lines[k])
    k += len(itf_files_fl)

# Faden num (first Faden num line)
if plot_Faden:
    legend_handles.append(lines[k])
    k += 1
    # Faden exp (first Faden exp line)
    legend_handles.append(lines[k])

plt.legend(legend_handles, ['Partitioned', 'Fluent', 'Faden num', 'Faden exp'], fontsize=14)

plt.tight_layout()
plt.savefig('interface_figures/interfaces_multiple_times.png', dpi=150)
plt.show()
plt.close()