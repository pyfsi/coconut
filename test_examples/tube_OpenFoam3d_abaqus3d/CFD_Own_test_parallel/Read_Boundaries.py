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

'''Piece of code to test whether face coordinates of the patches can be extracted simply from the boundary file in combination with the faces file'''

b_file = './constant/polyMesh/boundary'

f_file = './constant/polyMesh/faces'

interface = ["mantle"]


for name in interface:
    f_centres = []
    with open(b_file) as b:
        bool = 0
        for line in b:
            if name in line:
                nFaces=-1
                startFace=-1
                while(nFaces==-1 or startFace==-1):
                    l = b.readline()
                    if "nFaces" in l:
                        nFaces = int(l[:-2].split()[1])
                    if "startFace" in l:
                        startFace = int(l[:-2].split()[1])
    with open(f_file) as f:



