#Import of utilities
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.pyplot import xticks, yticks
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
import sys
import os  # to be able to run Linux terminal commands
from matplotlib import cm
import warnings
import pickle
import matplotlib.animation as animation

matplotlib.rcParams['font.size'] = 20
matplotlib.rcParams['axes.linewidth'] = 1
matplotlib.rcParams['lines.linewidth'] = 2
matplotlib.rcParams['axes.labelweight'] = 'bold'
matplotlib.rcParams['xtick.major.size'] = 6
matplotlib.rcParams['xtick.minor.size'] = 2.4
matplotlib.rcParams['ytick.major.size'] = 6
matplotlib.rcParams['ytick.minor.size'] = 2.4
matplotlib.rcParams['savefig.dpi'] = 300

def find_string_in_file(self, string_name, file_name):
    index = -1
    with open(file_name) as f:
        for num, line in enumerate(f):
            if string_name in line:
                index = num
                break
    f.close()
    return index

#### --------------- #####

directory = "./CFD_test5/processor0/0/"

X=[]
Y=[]
Z=[]

for name in ["ccx","ccy","ccz"]:
    file = directory+name
    with open(file) as f:
        line = f.readline()
        while not "mantle" in line:
            line=f.readline()
        for i in range(0,6):
            line=f.readline()
        while not ")" in line:
            if name=="ccx":
                X.append(float(line))
            elif name =="ccy":
                Y.append(float(line))
            elif name =="ccz":
                Z.append(float(line))
            line = f.readline()

R= []
for i in range(0,len(X)):
    R.append(np.sqrt(Y[i]**2+Z[i]**2))
R = np.array(R)
print(R.mean())
plt.plot(range(0,len(X)),R)

plt.show()