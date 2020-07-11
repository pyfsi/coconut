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

file = "./constant/polyMesh/points"

ORIG = []
with open(file) as orig:
    line = orig.readline()
    while not "(" in line:
        line=orig.readline()
    line=orig.readline()
    while any(char.isdigit() for char in line):
        coords = line[1:-2].split()
        ORIG.append([float(coords[0]),float(coords[1]),float(coords[2])])
        line=orig.readline()
ORIG=np.array(ORIG)
N=ORIG.shape[0]


file = "./1e-05/polyMesh/points"
NEW = []
with open(file) as new:
    line = new.readline()
    while not "(" in line:
        line=new.readline()
    line=new.readline()
    while any(char.isdigit() for char in line):
        coords = line[1:-2].split()
        NEW.append([float(coords[0]),float(coords[1]),float(coords[2])])
        line=new.readline()
NEW=np.array(NEW)


### nodes part of surface ###
fb = './constant/polyMesh/boundary'
ff = './constant/polyMesh/faces'

with open(fb) as f:
    line = f.readline()
    while not "mantle" in line:
        line = f.readline()
    for i in range(0, 4):
        line = f.readline()
    nFaces = int(line[:-2].split()[1])
    line = f.readline()
    startFace = int(line[:-2].split()[1])
print(nFaces)
print(startFace)

Fnodes = []
with open(ff) as f:
    line = f.readline()
    while not "(" in line:
        line = f.readline()
    line=f.readline()
    while any(char.isdigit() for char in line):
        list = line[2:-2].split()
        Fnodes.append(list)
        line = f.readline()

print(Fnodes[startFace])
print(Fnodes[startFace+1])
print(Fnodes[startFace+2])
print(Fnodes[startFace+3])

print(nFaces)
print(startFace)


plt.plot(range(0,N),NEW[:,1]-ORIG[:,1])

plt.show()