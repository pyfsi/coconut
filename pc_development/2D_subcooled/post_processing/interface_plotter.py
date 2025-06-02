import matplotlib.pyplot as plt
import numpy as np
import io
from scipy.signal import savgol_filter
from coconut.examples.post_processing.post_processing import *

# Start delay of simulation due to starting with LF = 0.01 & 0.02 (round to 0.1)
t_001 = 12.5 # s (based on conduction only Fluent simulation)
t_002 = 61.2 # s (based on CoCoNuT simulations starting from LF = 0.02)
t_restart = t_002 + 1800.0 # time of restart after manual remeshing at t = 1800 s

# different cases to be plotted
common_path = '../'
case_paths = ['Faden_split_3/case_results.pickle', 'Faden_split_3/run_1/case_results.pickle', 'Faden_split_3/HPC/case_results.pickle']
legend_entries = ['No initial sizing - remeshed', 'No initial sizing - not remeshed', 'No initial sizing - HPC']
dt = [0.1, 0.1, 0.1] # s
t_delay = [t_restart, t_002, t_002]

# fluent interfaces
plot_fluent = False
common_path_fl = './fluent_interfaces/'
itf_file_fl = 'itf-pos-1800-00s.xy'

parts = itf_file_fl.split('-')
legend_fl = 'Fluent - ' + parts[2] + '.' + parts[3].replace('s.xy', '') + ' s'
time_fl = float(parts[2]) + float(parts[3].replace('s.xy', ''))/100

# Paper Faden interfaces
plot_Faden = False
common_path_Fa = './Faden_paper/'
itf_file_Fa = 'Faden-num-itf-1800s.csv'

parts = itf_file_Fa.split('-')
legend_Fa = 'Faden - ' + parts[3].replace('s.csv', '') + ' s'

line_styles = ['r--', 'g--', 'b--', 'k--', 'r--', 'k--']

# Compare interface at a certain time: t_sim is the simulation time of the simulation with the largest time delay
t_sim = 90.0

if plot_fluent:
    if all(x == t_delay[0] for x in t_delay):
        t_sim = round(time_fl - t_delay[0], 1)
        print('Given time overruled by time of the Fluent interface.')
    else:
        print('Error: not all time offsets in t_delay are equal. Defaulting to given t_sim.')

t_delay_array = np.array(t_delay)
t_delta = np.max(t_delay_array) - t_delay_array
time = t_delta + t_sim

lines = []

for j, file in enumerate(case_paths):
    try:
        pp = PostProcess(common_path + file)
    except FileNotFoundError as e:
        print(f"Error: File not found - {e.filename}. File will be skipped.")
        continue

    sx = pp.add_subset(interface='interface_x', model_part='boundary_in_nodes')
    sy = pp.add_subset(interface='interface_y', model_part='boundary_out_faces')

    x = sx.get_values('coordinates', 'x')
    y = sx.get_values('coordinates', 'y')

    try:
        x = x[int(time[j] / dt[j]), :].flatten()
        y = y[int(time[j] / dt[j]), :].flatten()
        line, = plt.plot(x, y, line_styles[j], label=legend_entries[j])
        lines.append(line)

        print(f"Case '{legend_entries[j]}' plotted at total time {t_delay_array[j] + t_delta[j] + t_sim} s.")

    except IndexError:
        print(f"Case '{legend_entries[j]}' has not reached simulation time {time[j]} s yet (total time {np.max(t_delay_array) + t_sim} s).")

# fluent interface
if plot_fluent:
    with open(common_path_fl + itf_file_fl, 'r') as f:
        txt_lines = f.readlines()

    # Remove the first 4 header lines and the last line
    data_lines = txt_lines[4:-1]

    # Join the data lines back into a single string
    data_string = ''.join(data_lines)

    # Use numpy.loadtxt to efficiently read the data from the string.
    data = np.loadtxt(io.StringIO(data_string), delimiter='\t')
    y = data[:, 0]
    x = data[:, 1]
    line, = plt.plot(x, y, line_styles[len(case_paths)], label=legend_fl)
    lines.append(line)

# Faden interface
if plot_Faden:
    x, y = np.loadtxt(common_path_Fa + itf_file_Fa, skiprows=1, delimiter=',', unpack=True) # measured in mm
    sorted_indices = np.argsort(y)
    x = x[sorted_indices]
    y = y[sorted_indices]

    # --- Apply smoothing to the x data ---
    # Choose appropriate window_length and polyorder
    # window_length should be an odd integer, typically larger for more smoothing
    # polyorder should be less than window_length, e.g., 2 or 3 for most cases
    window_length = 11  # Example: Adjust as needed, must be odd
    polyorder = 2  # Example: Adjust as needed

    x_smoothed = savgol_filter(x, window_length, polyorder)

    line, = plt.plot(x_smoothed/1000, y/1000, line_styles[len(case_paths)+1], label=legend_Fa)
    lines.append(line)

plt.ylabel('y-coordinate [m]')
plt.xlabel('x-coordinate [m]')
plt.xlim((0, 0.04))
plt.legend(handles=lines)
plt.savefig('interface_figures/interface.png')
plt.show()
plt.close()