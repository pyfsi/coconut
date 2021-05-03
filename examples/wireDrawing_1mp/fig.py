
import matplotlib.pyplot as plt
import matplotlib
from coconut.coupling_components.component import Component
from coconut import data_structure

from subprocess import check_call
import numpy as np
import os
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

file_name = os.path.join('CSM', 'CSM_Time0Surface0Nodes.dat')
Nodes = np.genfromtxt (file_name)

file_name1 = os.path.join('CSM','CSM_Time1Surface0Output.dat' )
displacement = np.genfromtxt (file_name1)

file_name2 = os.path.join('CSM','CSM_Time1Surface0Output.dat')
displacementlast_iter=np.genfromtxt(file_name2)

file_name3 = os.path.join('CFD/postProcessing/PRESSURE_Wire/surface/0.010000', 'p_patch_Wire.raw')
pressure = np.genfromtxt (file_name3)

x = Nodes[:,1] + displacement[:,1]
x_axis = x[1:]
radial_axis = displacement[1:,0]

x1 = Nodes[:,1] + displacementlast_iter[:,1]
x1_axis = x1[1:]
radial1_axis = displacementlast_iter[1:,0]

fig, ax1 = plt.subplots()
# fig.suptitle('radial displacement wire drawing after 1 iteration', fontsize = 20)
fig.suptitle('Pressure distribution and radial displacement wire drawing after final iterations', fontsize = 20)

ax1.set_xlabel('wire position(m)')


color2 = 'tab:blue'
ax1.set_ylabel("Radial displacement surface wire (m)")
ax1.plot(x1_axis, radial1_axis, color=color2)
ax1.tick_params(axis='y', labelcolor=color2)

ax2 = ax1.twinx()

color = 'tab:red'
ax2.set_ylabel('Pressure fluid film(Pa)')
ax2.plot(pressure[:,0], pressure[:,3], color=color)
ax2.tick_params(axis = 'y', labelcolor=color)


fig.tight_layout()
plt.show()

