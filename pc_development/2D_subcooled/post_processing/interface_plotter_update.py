import matplotlib.pyplot as plt
import numpy as np
import io
from scipy.signal import savgol_filter
from coconut.examples.post_processing.post_processing import PostProcess

# ---------------- Settings ----------------
run = 'rerun'  # 'coarse' or 'rerun'

# Initial delay due to starting with LF < 0.02
t_001 = 12.5  # s
t_002 = 61.2  # s
t_run_1 = t_002  # initial offset

# Simulation times per run
if run == 'coarse':
    sim_times = [1800, 1620, 900, 1440, 720, 540, 180]  # coarse
elif run == 'rerun':
    sim_times = [1800, 1080, 1080, 1260, 1440, 5400]  # rerun

if run == 'coarse':
    common_path = '../Faden_split_3_coarse/'
elif run == 'rerun':
    common_path = '../Faden_split_3_rerun/'

# Time instances to plot
t_sims = [3600.0, 7200.0]  # seconds

# Partitioned simulation settings
dt = 0.1

# Fluent interface settings
plot_fluent = True
common_path_fl = './fluent_interfaces/'
itf_files_fl = ['itf-pos-3600-00s.xy', 'itf-pos-7200-00s.xy']

# Faden paper settings
plot_Faden = True
common_path_Fa = './Faden_paper/'
line_styles = ['r--', 'b--', 'r--', 'b--']  # different style exp./num.
itf_files_Fa = ['Faden-num-itf-3600s.csv', 'Faden-exp-itf-3600s.csv',
                'Faden-num-itf-7200s.csv', 'Faden-exp-itf-7200s.csv']

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
        line, = plt.plot(x_smooth/1000, y/1000, line_styles[i % len(line_styles)])
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

# Partitioned (take the first Partitioned line as representative)
legend_handles.append(lines[0])
# Fluent (first Fluent line)
legend_handles.append(lines[len(t_sims)])
# Faden num (first Faden num line)
legend_handles.append(lines[len(t_sims)*2])
# Faden exp (first Faden exp line)
legend_handles.append(lines[len(t_sims)*2 + 1])

plt.legend(legend_handles, ['Partitioned', 'Fluent', 'Faden num', 'Faden exp'], fontsize=14)

plt.tight_layout()
plt.savefig('interface_figures/interfaces_multiple_times.png', dpi=150)
plt.show()
plt.close()