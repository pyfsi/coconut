import matplotlib.pyplot as plt
import numpy as np
from coconut.examples.post_processing.post_processing import *

# different cases to be plotted
common_path = '../'
case_paths = ['Faden_full_5/case_results.pickle', 'Faden_full_6/case_results.pickle', 'Faden_full_7/case_results.pickle']
legend_entries = ['PRESTO!', 'bfw', 'bfw - p out']
dt = [0.1, 0.05, 0.05] # s

# fluent interfaces
plot_fluent = False
common_path_fl = './fluent_interfaces/'
itf_file = 'itf-738-97s.xy'
legend_fl = 'Fluent - 738.97 s'

line_styles = ['r--', 'g--', 'b--', 'k--', 'c--']

# Compare interface at a certain time
t = 30 # s

lines = []

for j, file in enumerate(case_paths):
    pp = PostProcess(common_path + file)

    sx = pp.add_subset(interface='interface_x', model_part='boundary_in_nodes')
    sy = pp.add_subset(interface='interface_y', model_part='boundary_out_faces')

    x = sx.get_values('coordinates', 'x')
    y = sx.get_values('coordinates', 'y')

    try:
        x = x[int(t / dt[j]), :].flatten()
        y = y[int(t / dt[j]), :].flatten()
        line, = plt.plot(x, y, line_styles[j], label=legend_entries[j])
        lines.append(line)
    except IndexError:
        print(f"Case '{legend_entries[j]}' has not reached time {t} s yet.")

# fluent interface
if plot_fluent:
    data = np.loadtxt(common_path_fl + itf_file, delimiter='\t')
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