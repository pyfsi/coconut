# Victor Van Riet
# 20/05/2024
#
# Python-program to plot nice figures from Fluent *.out-files
# Make sure the desired *.out file is present in this directory

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

debug = True
reverse = False

# different cases to be plotted
common_path = '/cfdfile2/data/fm/victor/Software/coconut/thermal_development/Eccomas/'
f = '/CFD_2/report-file.out'

if debug:
    case_paths = ['FFTB_relaxation' + f, 'TFFB_aitken' + f, 'FFTB_aitken' + f, 'TFFB_iqni' + f, 'FFTB_iqni' + f]
    legend_entries = ['FFTB_relaxation', 'TFFB_aitken', 'FFTB_aitken', 'TFFB_iqni', 'FFTB_iqni']
else:
    case_paths = ['TFFB_relaxation' + f, 'FFTB_relaxation' + f, 'TFFB_aitken' + f, 'FFTB_aitken' + f, 'TFFB_iqni' + f, 'FFTB_iqni' + f]
    legend_entries = ['TFFB_relaxation', 'FFTB_relaxation', 'TFFB_aitken', 'FFTB_aitken', 'TFFB_iqni', 'FFTB_iqni']

line_styles = ['b--', 'b.-', 'g--', 'g.-', 'r--', 'r.-']
line_temp = []
line_hf = []
variables = ['temperature', 'heat-flux']
units = ['[K]', '[W/m^2]']

for i, name in enumerate(legend_entries):
    # Collect only the numerical data first
    if reverse:
        full_path = common_path + 'Transient/reverse/' + case_paths[i]
    else:
        full_path = common_path + 'Transient/' + case_paths[i]
    csv_path = common_path + 'post_processing/csv/' + name + '.csv'
    read_file = pd.read_csv(full_path, delimiter='\s+', skiprows=[0, 1, 2])
    try:
        os.remove(csv_path)
    except:
        print('No *.csv file present. Check if the file is closed!')

for j, var in enumerate(variables):
    for i, name in enumerate(legend_entries):
        csv_path = common_path + 'post_processing/csv/' + name + '.csv'
        read_file.to_csv(csv_path, index=None)

        # Create an array of the data csv
        data_array = np.loadtxt(csv_path, delimiter=',')

        time = data_array[0:-1,3]
        if var == 'temperature':
            data = data_array[0:-1,2]
        else:
            data = data_array[0:-1,1]

        # Create lines
        line, = plt.plot(time, data, line_styles[i], label=name)
        if var == 'temperature':
            line_temp.append(line)
        else:
            line_hf.append(line)

    # Plot lines
    plt.ylabel('Avg. interface ' + var + ' ' + units[j])
    plt.xlabel('Time [s]')
    plt.legend(handles=line_temp)
    plt.savefig('./figures/avg-itf-' + var + '-transient.svg')
    plt.show()
    plt.close()