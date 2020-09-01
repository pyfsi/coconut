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


#parallel setting
nproc = 2

def Get_Transform(dir,orig_boundary_Ind):
    f_p = f"./{dir}/constant/polyMesh/points"  # To read the number of points
    f_b = f'./{dir}/constant/polyMesh/boundary'
    f_f = f'./{dir}/constant/polyMesh/faces'

    ### For sequence of pointDisplacement file ###
    # Identify face numbers that belong to the desired boundary
    name = 'mantle'
    with open(f_b) as f:
        line = f.readline()
        while not name in line:
            line = f.readline()
        for i in range(0, 4):
            line = f.readline()
        nFaces = int(line[:-2].split()[1])
        line = f.readline()
        startFace = int(line[:-2].split()[1])
    print(nFaces)
    print(startFace)

    # Identify points (in sequence of pointDisplacement file) that belong to the boundary
    Fnodes = []
    with open(f_f) as f:
        line = f.readline()
        while not "(" in line:
            line = f.readline()
        line = f.readline()
        while any(char.isdigit() for char in line):
            list = line[2:-2].split()
            Fnodes.append(list)
            line = f.readline()

    # Get number of points to keep a list of booleans
    for _ in (True,):  # Creates a loop which is executed once, but can be exited with a break statement
        prev_line = "("
        with open(f_p) as f:
            line = f.readline()
            while not "(" in line:
                prev_line = line
                line = f.readline()
            break
    N_points = int(prev_line)
    print(N_points)

    points_Bool = np.zeros((N_points, 1))  # !! Carefull indices in OpenFoam start at 1
    boundary_Ind = []
    All_Fnodes = []
    with open(f_f) as f:
        line = f.readline()
        while not "(" in line:
            line = f.readline()
        line = f.readline()
        while any(char.isdigit() for char in line):
            list = line[2:-2].split()
            All_Fnodes.append(list)
            line = f.readline()

    for nf in range(startFace, startFace + nFaces):
        list = All_Fnodes[nf]
        for n in list:
            if points_Bool[int(n) - 1, 0] == 0:  # If not yet in list, add it to the list
                boundary_Ind.append(int(n))
                points_Bool[int(n) - 1, 0] = 1
    boundary_Ind = np.array(boundary_Ind)

    f_swap = f"./{dir}/constant/polyMesh/pointProcAddressing" #need this to associate the correct indices
    seq = []
    with open(f_swap) as f:
        line = f.readline()
        while not "(" in line:
            line = f.readline()
        line = f.readline()
        while any(char.isdigit() for char in line):
            ind = int(line)
            seq.append(ind)
            line = f.readline()

    swap = []
    for ind in boundary_Ind:
        x = seq[ind]
        X = np.where(orig_boundary_Ind == x+1)
        print(X)
        swap.append(X[0][0])

    test = 1

    return swap

##

'''Write a pointDisplacement file to test the displacement of the points based on their coordinates in parallel
Can the updated pointDisplacement for later time steps only list those for the mantle?'''

f_p = "./constant/polyMesh/points" #To read the number of points
f_b = './constant/polyMesh/boundary'
f_f = './constant/polyMesh/faces'

### For sequence of pointDisplacement file ###
# Identify face numbers that belong to the desired boundary
name = 'mantle'
with open(f_b) as f:
    line = f.readline()
    while not name in line:
        line = f.readline()
    for i in range(0, 4):
        line = f.readline()
    nFaces = int(line[:-2].split()[1])
    line = f.readline()
    startFace = int(line[:-2].split()[1])
print(nFaces)
print(startFace)

# Identify points (in sequence of pointDisplacement file) that belong to the boundary
Fnodes = []
with open(f_f) as f:
    line = f.readline()
    while not "(" in line:
        line = f.readline()
    line=f.readline()
    while any(char.isdigit() for char in line):
        list = line[2:-2].split()
        Fnodes.append(list)
        line = f.readline()

#Get number of points to keep a list of booleans
for _ in (True,): #Creates a loop which is executed once, but can be exited with a break statement
    prev_line = "("
    with open(f_p) as f:
        line = f.readline()
        while not "(" in line:
            prev_line=line
            line = f.readline()
        break
N_points=int(prev_line)
print(N_points)

points_Bool = np.zeros((N_points,1)) #!! Carefull indices in OpenFoam start at 1
boundary_Ind = []
All_Fnodes = []
with open(f_f) as f:
    line = f.readline()
    while not "(" in line:
        line = f.readline()
    line=f.readline()
    while any(char.isdigit() for char in line):
        list = line[2:-2].split()
        All_Fnodes.append(list)
        line = f.readline()

for nf in range(startFace,startFace+nFaces):
        list = All_Fnodes[nf]
        for n in list:
            if points_Bool[int(n),0]==0: #If not yet in list, add it to the list
                boundary_Ind.append(int(n))
                points_Bool[int(n),0]=1
boundary_Ind = np.array(boundary_Ind)

## For original coordinates of the points ###
n_points_b = len(boundary_Ind)

print(f"Number of mesh point on boundary {name}:{n_points_b}")
ORIG = []
    # get points and their coordinates
with open(f_p) as f:
    line = f.readline()
    while not "(" in line:
        line = f.readline()
    line = f.readline()
    while any(char.isdigit() for char in line):
        coords = line[1:-2].split()
        ORIG.append([float(coords[0]), float(coords[1]), float(coords[2])])
        line = f.readline()
ORIG = np.array(ORIG)
orig_coords_b=[]
for b in boundary_Ind:
    # print(b)
    orig_coords_b.append(ORIG[b,:])
orig_coords_b = np.array(orig_coords_b)

print("Launch TS 1")
os.system("sed -i 's/.*startTime.*/startTime 0;/' ./system/controlDict  ") #Start from 0
os.system("sed -i 's/.*endTime.*/endTime 0.00001;/' ./system/controlDict  ")
## Impose a sinusoidal variation of the radius for the tube
R_orig = 0.005
L = 0.05

delta_R = R_orig/10


disp_b = np.zeros(np.shape(orig_coords_b))
for i in range(0,orig_coords_b.shape[0]):
    r_vec = np.array([0,orig_coords_b[i,1],orig_coords_b[i,2]])
    r_vec_n = r_vec/np.linalg.norm(r_vec)
    X = orig_coords_b[i,0]
    dX = 0
    dY = delta_R*np.sin(2*np.pi*(X+L/2.0)/L)*r_vec_n[1]
    dZ = delta_R*np.sin(2*np.pi*(X+L/2.0)/L)*r_vec_n[2]
    disp_b[i,0] = dX
    disp_b[i, 1] = dY
    disp_b[i, 2] = dZ

new_coords_b=orig_coords_b+disp_b

plt.figure()
plt.scatter(new_coords_b[:,0],new_coords_b[:,1])

plt.figure()
plt.scatter(new_coords_b[:,0],new_coords_b[:,2])

## Create first pointdisplacement file ##
os.system("cp -r /lusers/temp/lucas/Temp_CoCoNut/coconut/coupling_components/solver_wrappers/openfoam/pointDisplacement_raw ./0/pointDisplacement")
f_disp = "./0/pointDisplacement"

boundaries = [name]

for bound in boundaries:
    with open(f_disp,"a") as f:
        f.write(bound+"\n")
        f.write("{\n")
        f.write("type fixedValue;\n")
        f.write("value  nonuniform List<vector>\n")
        f.write(f"{n_points_b}\n")
        f.write(f"(\n")
        for i in range(0,disp_b.shape[0]):
            f.write(f"({disp_b[i,0]} {disp_b[i,1]} {disp_b[i,2]})\n")
        f.write(");\n")
        f.write("}\n")

os.system("decomposePar")
os.system(f"mpirun -np {nproc} pimpleDyMFoam -parallel &> log_first_TS")


## second time step in parallel
disp_b = np.zeros(np.shape(orig_coords_b)) #original serial sequence
for i in range(0,orig_coords_b.shape[0]):
    r_vec = np.array([0,orig_coords_b[i,1],orig_coords_b[i,2]])
    r_vec_n = r_vec/np.linalg.norm(r_vec)
    X = orig_coords_b[i,0]
    dX = 0
    dY = delta_R*np.sin(2*np.pi*(X+L/2.0)/L)*r_vec_n[1]
    dZ = delta_R*np.sin(2*np.pi*(X+L/2.0)/L)*r_vec_n[2]
    disp_b[i,0] = dX
    disp_b[i, 1] = dY
    disp_b[i, 2] = dZ

new_coords_b=orig_coords_b+disp_b

print("Start TS 2")
os.system("sed -i 's/.*startTime.*/startTime 0.00001;/' ./system/controlDict ") #Start from 0
os.system("sed -i 's/.*endTime.*/endTime 0.00002;/' ./system/controlDict  ")

# proc_directories = [f"processor{i}" for i in range(0,nproc)]
#
# for dir in proc_directories:
#     Get_Transform(dir,boundary_Ind)


os.system("reconstructPar -time '1e-05' -fields '(pointDisplacement)'")

os.system("mv ./1e-05/pointDisplacement  ./1e-05/pointDisplacement_orig")

delta_R *= -1
#reverse the deformation
disp_b = np.zeros(np.shape(orig_coords_b))
for i in range(0,orig_coords_b.shape[0]):
    r_vec = np.array([0,orig_coords_b[i,1],orig_coords_b[i,2]])
    r_vec_n = r_vec/np.linalg.norm(r_vec)
    X = orig_coords_b[i,0]
    dX = 0
    dY = delta_R*np.sin(2*np.pi*(X+L/2.0)/L)*r_vec_n[1]
    dZ = delta_R*np.sin(2*np.pi*(X+L/2.0)/L)*r_vec_n[2]
    disp_b[i,0] = dX
    disp_b[i, 1] = dY
    disp_b[i, 2] = dZ

new_coords_b=orig_coords_b+disp_b

f_disp = "./1e-05/pointDisplacement"

os.system(f"cp -r /lusers/temp/lucas/Temp_CoCoNut/coconut/coupling_components/solver_wrappers/openfoam/pointDisplacement_raw {f_disp}")

boundaries = [name]

for bound in boundaries:
    with open(f_disp,"a") as f:
        f.write(bound+"\n")
        f.write("{\n")
        f.write("type fixedValue;\n")
        f.write("value  nonuniform List<vector>\n")
        f.write(f"{n_points_b}\n")
        f.write(f"(\n")
        for i in range(0,disp_b.shape[0]):
            f.write(f"({disp_b[i,0]} {disp_b[i,1]} {disp_b[i,2]})\n")
        f.write(");\n")
        f.write("}\n")

os.system("decomposePar -fields -time '1e-05'")
os.system(f"mpirun -np {nproc} pimpleDyMFoam -parallel &> log_second_TS")
# os.system("pimpleDyMFoam &> log_second_TS")
#
# plt.figure(1)
# plt.scatter(new_coords_b[:,0],new_coords_b[:,1])
#
# plt.figure(2)
# plt.scatter(new_coords_b[:,0],new_coords_b[:,2])
#
#
# ## Do some more iterations to check if quality remained ok
# os.system("sed -i 's/.*startTime.*/startTime 0.00002;/' ./system/controlDict ") #Start from 0
# os.system("sed -i 's/.*endTime.*/endTime 0.0001;/' ./system/controlDict  ")
#
# os.system("pimpleDyMFoam &> log_extra_TS")


plt.show()