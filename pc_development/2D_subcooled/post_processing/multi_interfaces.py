import matplotlib.pyplot as plt
import numpy as np
import io
from scipy.signal import savgol_filter
from coconut.examples.post_processing.post_processing import *

# Start delay of simulation due to starting with LF = 0.01 & 0.02 (round to 0.1)
t_001 = 12.5 # s (based on conduction only Fluent simulation)
t_002 = 61.2 # s (based on CoCoNuT simulations starting from LF = 0.01)

# Coconut case to compare with Fluent interfaces
common_path = '../'
case_path = 'Faden_split_3/case_results.pickle'
dt = 0.1 # s
t_delay = t_002

pp = PostProcess(common_path + case_path)
sx = pp.add_subset(interface='interface_x', model_part='boundary_in_nodes')
sy = pp.add_subset(interface='interface_y', model_part='boundary_out_faces')

x_coco = sx.get_values('coordinates', 'x')
y_coco = sx.get_values('coordinates', 'y')

# Fluent interfaces
common_path_fl = './fluent_interfaces/'
itf_file_fl = ['itf-pos-257-96s.xy', 'itf-pos-738-97s.xy', 'itf-pos-1200-00s.xy', 'itf-pos-1800-00s.xy']
legend = []
times_fl_rounded = []
for itf in itf_file_fl:
    parts = itf.split('-')
    legend.append(parts[2] + '.' + parts[3].replace('s.xy', '') + ' s')
    times_fl_rounded.append(round(float(parts[2]) + float(parts[3].replace('s.xy', '')) / 100, 1))

# To be filled with each plot line
lines = []

for j, t in enumerate(times_fl_rounded):
    try:
        x = x_coco[int((t - t_delay) / dt), :].flatten()
        y = y_coco[int((t - t_delay) / dt), :].flatten()
        line, = plt.plot(x, y, 'k-', label='Partitioned - ' + legend[j])
        lines.append(line)

    except IndexError:
        print(f"Coconut has not reached simulation time {t} s yet.")

    # Fluent interface
    with open(common_path_fl + itf_file_fl[j], 'r') as f:
        txt_lines = f.readlines()

    # Remove the first 4 header lines and the last line
    data_lines = txt_lines[4:-1]

    # Join the data lines back into a single string
    data_string = ''.join(data_lines)

    # Use numpy.loadtxt to efficiently read the data from the string.
    data = np.loadtxt(io.StringIO(data_string), delimiter='\t')
    y = data[:, 0]
    x = data[:, 1]
    line, = plt.plot(x, y, 'r--', label='Enthalpy-Porosity - ' + legend[j])
    lines.append(line)

# Paper Faden interfaces
plot_Faden = True
common_path_Fa = './Faden_paper/'
itf_file_Fa = 'Faden-num-itf-1800s.csv'

parts = itf_file_Fa.split('-')
legend_Fa = 'Faden - ' + parts[3].replace('s.csv', '') + ' s'

# Faden interface
if plot_Faden:
    x, y = np.loadtxt(common_path_Fa + itf_file_Fa, skiprows=1, delimiter=',', unpack=True)  # measured in mm
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

    line, = plt.plot(x_smoothed / 1000, y / 1000, 'g--', label=legend_Fa)
    lines.append(line)

# Also compare Fluent & Faden at 1800 s & 3600 s (both num. and exp.) in the same plot to give an idea about accuracy of own Fluent sim.

plt.ylabel('y-coordinate [m]', fontsize=16)
plt.xlabel('x-coordinate [m]', fontsize=16)
plt.xlim((0, 0.014))
plt.xticks(fontsize=14)  # Increase font size of x-ticks
plt.yticks(fontsize=14)  # Increase font size of y-ticks
#plt.legend(handles=lines, fontsize=16)
plt.tight_layout()
plt.savefig('interface_figures/multi_interfaces.png', dpi=150)
plt.show()
plt.close()