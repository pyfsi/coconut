import matplotlib.pyplot as plt
import matplotlib
from coconut.coupling_components.component import Component
from coconut import data_structure

from subprocess import check_call
import numpy as np
import os
from scipy import interpolate
import sys
import time
import subprocess
import re

matplotlib.rcParams['font.size'] = 14
matplotlib.rcParams['axes.linewidth']=1
matplotlib.rcParams['lines.linewidth'] = 2
matplotlib.rcParams['axes.labelweight'] = 'bold'
matplotlib.rcParams['xtick.major.size'] = 4
matplotlib.rcParams['xtick.minor.size'] = 2.4
matplotlib.rcParams['ytick.major.size'] = 4
matplotlib.rcParams['ytick.minor.size'] = 2.4
matplotlib.rcParams['savefig.dpi'] = 300

X_name = os.path.join('CFD', '0/Cx')
Y_name = os.path.join('CFD', '0/Cy')
Z_name = os.path.join('CFD', '0/Cz')


timestamp = '0.00330'
number_of_iterations = 3

pointDisplacement = os.path.join('CFD',timestamp, 'cellDisplacement')


pointDisplacement_fluid = []
pointDisplacement_structure = []

for i in range(number_of_iterations):
    file_name_displacement_structure = os.path.join('CSM/postProcessing/POINT_DISP_wireTopContact/surface/',timestamp,
                                                    f'dispPoint_{i+1}/dispPoint')
    with open(file_name_displacement_structure) as f:
        lines = f.readlines()
    array_plot = int(len(lines)/2)
    array_disp_struc = np.zeros((array_plot,3))
    k = 0
    for line in range(len(lines)):
        if (line % 2) == 0:
            coor = lines[line]
            for index in range(len(coor)):
                    if ' ' in coor[index]:
                        if index < 15:
                            array_disp_struc[k,0] = coor[0:index]
                            tmp = index
                        else:
                            array_disp_struc[k,1] = coor[tmp:index]
                            tmp2 = index
                    if ')' in coor[index]:
                        array_disp_struc[k,2] = coor[tmp2:index]
            k+=1
    # array_disp_struc[:,1] = array_disp_struc[:,1]*1000000
    pointDisplacement_structure.append(array_disp_struc)
    
    if i == 0:
        file_name_displacement_fluid = os.path.join('CFD',timestamp,f'pointDisplacement')
    else:
        file_name_displacement_fluid = os.path.join('CFD', timestamp,f'pointDisplacement_Next_Iter{i}')
    x_axis_file = os.path.join('CFD/postProcessing/mp_x0_points')
    fTemp = open(file_name_displacement_fluid, 'r')
    fTempLines = fTemp.readlines()
    for line in range(len(fTempLines)):
        if "Wire" in fTempLines[line]:
            count = int(fTempLines[line + 4])
            array_disp_fluid = np.zeros((count,3))
            for i in range(count):
                coor = (fTempLines[line + 5 + i])
                for index in range(len(coor)):
                    if ' ' in coor[index]:
                        if index < 25:
                            # print(coor[2:index])
                            array_disp_fluid[i,0] = 0.0
                            tmp = index
                        else:
                            array_disp_fluid[i,1] = coor[tmp:index]
                            tmp2 = index
    fPoint = open(x_axis_file,'r')
    fPointTemp = fPoint.readlines()
    for line in range(len(fPointTemp)):
      coor = fPointTemp[line]
      array_disp_fluid[line,0] = coor[0:len(coor)]

    array_plot_fluid = int(len(fPointTemp) / 2)

    array_disp_fluid_final = np.zeros((array_plot_fluid, 3))

    m = 0
    for line in range(array_disp_fluid[:,0].size):
        # print(array_disp[line,0])
        if (line % 2) == 0:
            array_disp_fluid_final[m,0] = array_disp_fluid[line,0]
            array_disp_fluid_final[m,1] = array_disp_fluid[line,1]
            m+=1
    pointDisplacement_fluid.append(array_disp_fluid_final)

fig1, ax1 = plt.subplots()

plot_last_struc = pointDisplacement_structure[-1]
plot_fluid = pointDisplacement_fluid[-1]

f = interpolate.interp1d(plot_fluid[:,0],plot_fluid[:,1],fill_value = 'extrapolate')
corr = f(plot_last_struc[:,0])
print('corr')
print(corr)
print("plot_fluid)")
print(plot_last_struc[:,1])
num_fault = abs(plot_last_struc[:,1] - corr)

ax1.plot(plot_last_struc[:, 0], num_fault, label = f'final_displacement_fluid')
ax1.legend()
fig1.suptitle(f'Transfer radial displacement timestamp: {timestamp}', fontsize=20)
ax1.set_ylabel("radial displacement (m)")
ax1.set_xlabel('wire position(m)')
ax1.tick_params(axis='y')
plt.grid()

plt.show()


