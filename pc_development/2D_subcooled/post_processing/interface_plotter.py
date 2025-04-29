import matplotlib.pyplot as plt
import numpy as np
import io
from coconut.examples.post_processing.post_processing import *

# Start delay of simulation due to starting with LF = 0.01 & 0.02 (round to 0.1)
t_001 = 12.5 # s (based on conduction only Fluent simulation)
t_002 = 61.2 # s (based on CoCoNuT simulations starting from LF = 0.01)

# different cases to be plotted
common_path = '../'
case_paths = ['Faden_full_8/case_results.pickle', 'Faden_full_11/case_results.pickle','Faden_full_12/case_results.pickle']
legend_entries = ['var. density', 'var. density - double', 'cst. density - corrected']
dt = [0.05, 0.05, 0.05, 0.05] # s
t_delay = [t_001, t_002, t_001, t_001]

# fluent interfaces
plot_fluent = False
common_path_fl = './fluent_interfaces/'
itf_file = 'itf-pos-1200-00s.xy'

parts = itf_file.split('-')
legend_fl = 'Fluent - ' + parts[2] + '.' + parts[3].replace('s.xy', '') + ' s'

line_styles = ['r--', 'g--', 'b--', 'k--', 'c--']

# Compare interface at a certain time: t_sim is the simulation time of the simulation with the largest time delay
t_sim = 230 # s

t_delay_array = np.array(t_delay)
t_delta = np.max(t_delay_array) - t_delay_array
time = t_delta + t_sim

lines = []

for j, file in enumerate(case_paths):
    pp = PostProcess(common_path + file)

    sx = pp.add_subset(interface='interface_x', model_part='boundary_in_nodes')
    sy = pp.add_subset(interface='interface_y', model_part='boundary_out_faces')

    x = sx.get_values('coordinates', 'x')
    y = sx.get_values('coordinates', 'y')

    try:
        x = x[int(time[j] / dt[j]), :].flatten()
        y = y[int(time[j] / dt[j]), :].flatten()
        line, = plt.plot(x, y, line_styles[j], label=legend_entries[j])
        lines.append(line)

        print(f"Case '{legend_entries[j]}' plotted at total time {t_delay_array[j] + t_delta[j] + t_sim} s yet.")

    except IndexError:
        print(f"Case '{legend_entries[j]}' has not reached simulation time {time[j]} s yet (total time {np.max(t_delay_array) + t_sim} s).")

# fluent interface
if plot_fluent:
    with open(common_path_fl + itf_file, 'r') as f:
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

plt.ylabel('y-coordinate [m]')
plt.xlabel('x-coordinate [m]')
plt.xlim((0, 0.04))
plt.legend(handles=lines)
#plt.savefig('interface.png')
plt.show()
plt.close()