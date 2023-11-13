# This is a simple driver for the RVE FFT solver
# the material should start from the virgin state,
# i.e. F = identity matrix

import mesostress_v2

import pickle

import os
os.environ["MKL_NUM_THREADS"] = "1" 
os.environ["NUMEXPR_NUM_THREADS"] = "1" 
os.environ["OMP_NUM_THREADS"] = "1" 

import numpy as np
import scipy.sparse.linalg as sp
import itertools

import sys
import time

import argparse

parser = argparse.ArgumentParser(description='Process some integers.')

parser.add_argument('--nCells', metavar='nCells', type=int, nargs='+',
                    help='The number of cells in each dimension, must be an odd number',required=True)
parser.add_argument('--avgFSeries', metavar='avgFSeries', nargs='+',
                    help='File location for the average deformation gradient tensor interpolation table',required=True)
parser.add_argument('--nIncrements', metavar='nIncrements', nargs='+',
                    help='The number of increments to take',required=True)
parser.add_argument('--relaxationFactor', metavar='relaxationFactor', nargs='+',
                    help='The relaxation factor to apply at each increment iteration',default=0.9)
parser.add_argument('--maxIts', metavar='maxIts', nargs='+',
                    help='The maximum number of Newton iterations',default=200)
parser.add_argument('--outputFile', metavar='outputFile', nargs='+',
                    help='',required=True)
parser.add_argument('--usingNumericalGradient', metavar='usingNumericalGradient', nargs='+',
                    help='Switch to use a central difference approximation to the tangent matrix',default=False)
parser.add_argument('--inputPkl', metavar='inputPkl', nargs='+',
                    help='If the state variables are to be initialised from a saved pkl file',required=False, default=None)

# not yet implemented
#parser.add_argument('--captureYieldStrengthEvolution', metavar='captureYieldStrengthEvolution', nargs='+',
 #                   help='This indicates whether to evaluate various yield strengths',default=False)



# parse the arguments
args = parser.parse_args()
nCells = int(args.nCells[0])
nIncrements = int(args.nIncrements[0])
avgFSeries_loc = args.avgFSeries[0]
outputFile_loc = args.outputFile[0]
try:
    usingNumericalGradient = args.usingNumericalGradient[0]
except:
    usingNumericalGradient = args.usingNumericalGradient
try:
    defaultRelaxationFactor = float(args.relaxationFactor[0])
except:
    defaultRelaxationFactor = args.relaxationFactor
try:
    maxIts = int(args.maxIts[0])
except:
    maxIts = args.maxIts

if nCells % 2 == 0:
    raise Exception('nCells must be an odd integer value')

if usingNumericalGradient:
    print('Using central difference approximation for the tangent matrix')


try:
    inputPkl = args.inputPkl[0]
    print('Reading initial state variables from file', inputPkl)    
except:
    pass



    
# setup the deformation gradient interpolation
t_series = []
F11_series = []
F12_series = []
F13_series = []
F21_series = []
F22_series = []
F23_series = []
F31_series = []
F32_series = []
F33_series = []

with open(avgFSeries_loc) as f:
    lines = f.readlines()
    for l in lines[1:-1]:
        l = l.replace('(',' ')
        l = l.replace(')',' ')
        split_line = np.asarray([float(x) for x in l.split()])

        t_series.append(split_line[0])
        F11_series.append(split_line[1])
        F12_series.append(split_line[2])
        F13_series.append(split_line[3])        
        F21_series.append(split_line[4])
        F22_series.append(split_line[5])
        F23_series.append(split_line[6])        
        F31_series.append(split_line[7])
        F32_series.append(split_line[8])
        F33_series.append(split_line[9])        
        
F11 = lambda t: np.interp(t, t_series, F11_series)
F12 = lambda t: np.interp(t, t_series, F12_series)
F13 = lambda t: np.interp(t, t_series, F13_series)
F21 = lambda t: np.interp(t, t_series, F21_series)
F22 = lambda t: np.interp(t, t_series, F22_series)
F23 = lambda t: np.interp(t, t_series, F23_series)
F31 = lambda t: np.interp(t, t_series, F31_series)
F32 = lambda t: np.interp(t, t_series, F32_series)
F33 = lambda t: np.interp(t, t_series, F33_series)

T_final = t_series[-1]


# # turn of warning for zero division
# # (which occurs in the linearization of the logarithmic strain)
np.seterr(divide='ignore', invalid='ignore')

# # ----------------------------------- GRID ------------------------------------

Nx     = nCells          # number of voxels in x-direction
Ny     = nCells         # number of voxels in y-direction
Nz     = nCells          # number of voxels in z-direction
shape  = [Nx,Ny,Nz]  # number of voxels as list: [Nx,Ny,Nz]

# ----------------------------- TENSOR OPERATIONS -----------------------------

# tensor operations / products: np.einsum enables index notation, avoiding loops
# e.g. ddot42 performs $C_ij = A_ijkl B_lk$ for the entire grid
trans2 = lambda A2   : np.einsum('ijxyz          ->jixyz  ',A2   )
ddot22 = lambda A2,B2: np.einsum('ijxyz  ,jixyz  ->xyz    ',A2,B2)
ddot42 = lambda A4,B2: np.einsum('ijklxyz,lkxyz  ->ijxyz  ',A4,B2)
ddot44 = lambda A4,B4: np.einsum('ijklxyz,lkmnxyz->ijmnxyz',A4,B4)
dot11  = lambda A1,B1: np.einsum('ixyz   ,ixyz   ->xyz    ',A1,B1)
dot22  = lambda A2,B2: np.einsum('ijxyz  ,jkxyz  ->ikxyz  ',A2,B2)
dot24  = lambda A2,B4: np.einsum('ijxyz  ,jkmnxyz->ikmnxyz',A2,B4)
dot42  = lambda A4,B2: np.einsum('ijklxyz,lmxyz  ->ijkmxyz',A4,B2)
dyad22 = lambda A2,B2: np.einsum('ijxyz  ,klxyz  ->ijklxyz',A2,B2)
dyad11 = lambda A1,B1: np.einsum('ixyz   ,jxyz   ->ijxyz  ',A1,B1)

# eigenvalue decomposition of 2nd-order tensor: return in convention i,j,x,y,z
# NB requires to swap default order of NumPy (in in/output)
def eig2(A2):
    swap1i    = lambda A1: np.einsum('xyzi ->ixyz ',A1)
    swap2     = lambda A2: np.einsum('ijxyz->xyzij',A2)
    swap2i    = lambda A2: np.einsum('xyzij->ijxyz',A2)
    vals,vecs = np.linalg.eig(swap2(A2))
    vals      = swap1i(vals)
    vecs      = swap2i(vecs)
    return vals,vecs

# logarithm of grid of 2nd-order tensors
def ln2(A2):
    vals,vecs = eig2(A2)
    return sum([np.log(vals[i])*dyad11(vecs[:,i],vecs[:,i]) for i in range(3)])

# exponent of grid of 2nd-order tensors
def exp2(A2):
    vals,vecs = eig2(A2)
    return sum([np.exp(vals[i])*dyad11(vecs[:,i],vecs[:,i]) for i in range(3)])

# determinant of grid of 2nd-order tensors
def det2(A2):
    return (A2[0,0]*A2[1,1]*A2[2,2]+A2[0,1]*A2[1,2]*A2[2,0]+A2[0,2]*A2[1,0]*A2[2,1])-           (A2[0,2]*A2[1,1]*A2[2,0]+A2[0,1]*A2[1,0]*A2[2,2]+A2[0,0]*A2[1,2]*A2[2,1])

# inverse of grid of 2nd-order tensors
def inv2(A2):
    A2det = det2(A2)
    A2inv = np.empty([3,3,Nx,Ny,Nz])
    A2inv[0,0] = (A2[1,1]*A2[2,2]-A2[1,2]*A2[2,1])/A2det
    A2inv[0,1] = (A2[0,2]*A2[2,1]-A2[0,1]*A2[2,2])/A2det
    A2inv[0,2] = (A2[0,1]*A2[1,2]-A2[0,2]*A2[1,1])/A2det
    A2inv[1,0] = (A2[1,2]*A2[2,0]-A2[1,0]*A2[2,2])/A2det
    A2inv[1,1] = (A2[0,0]*A2[2,2]-A2[0,2]*A2[2,0])/A2det
    A2inv[1,2] = (A2[0,2]*A2[1,0]-A2[0,0]*A2[1,2])/A2det
    A2inv[2,0] = (A2[1,0]*A2[2,1]-A2[1,1]*A2[2,0])/A2det
    A2inv[2,1] = (A2[0,1]*A2[2,0]-A2[0,0]*A2[2,1])/A2det
    A2inv[2,2] = (A2[0,0]*A2[1,1]-A2[0,1]*A2[1,0])/A2det
    return A2inv


# ------------------------ INITIATE (IDENTITY) TENSORS ------------------------

# identity tensor (single tensor)
i    = np.eye(3)
# identity tensors (grid)
I    = np.einsum('ij,xyz'           ,                  i   ,np.ones([Nx,Ny,Nz]))
I4   = np.einsum('ijkl,xyz->ijklxyz',np.einsum('il,jk',i,i),np.ones([Nx,Ny,Nz]))
I4rt = np.einsum('ijkl,xyz->ijklxyz',np.einsum('ik,jl',i,i),np.ones([Nx,Ny,Nz]))
I4s  = (I4+I4rt)/2.
II   = dyad22(I,I)


# ------------------------------------ FFT ------------------------------------

# projection operator (only for non-zero frequency, associated with the mean)
# NB: vectorized version of "hyper-elasticity.py"
# - allocate / support function
Ghat4  = np.zeros([3,3,3,3,Nx,Ny,Nz])                # projection operator
x      = np.zeros([3      ,Nx,Ny,Nz],dtype='int64')  # position vectors
q      = np.zeros([3      ,Nx,Ny,Nz],dtype='int64')  # frequency vectors
delta  = lambda i,j: np.float(i==j)                  # Dirac delta function
# - set "x" as position vector of all grid-points   [grid of vector-components]
x[0],x[1],x[2] = np.mgrid[:Nx,:Ny,:Nz]
# - convert positions "x" to frequencies "q"        [grid of vector-components]
for i in range(3):
    freq = np.arange(-(shape[i]-1)/2,+(shape[i]+1)/2,dtype='int64')
    q[i] = freq[x[i]]
# - compute "Q = ||q||", and "norm = 1/Q" being zero for the mean (Q==0)
#   NB: avoid zero division
q       = q.astype(np.float)
Q       = dot11(q,q)
Z       = Q==0
Q[Z]    = 1.
norm    = 1./Q
norm[Z] = 0.
# - set projection operator                                   [grid of tensors]
for i, j, l, m in itertools.product(range(3), repeat=4):
    Ghat4[i,j,l,m] = norm*delta(i,m)*q[j]*q[l]

# # STRESS CONTROL
# for x in range(0,Nx):
#     for y in range(0,Ny):
#         for z in range(0,Nz):
#             for k in range(0,3):
#                 for l in range(0,3):
#                     if q[0,x,y,z] == q[1,x,y,z] == q[2,x,y,z] == 0:
#                         Ghat4[1,1,k,l] = delta(1,k)*delta(1,l)
#                         Ghat4[2,2,k,l] = delta(2,k)*delta(2,l)                        
# # END STRESS CONTROL
    
# (inverse) Fourier transform (for each tensor component in each direction)
fft  = lambda x: np.fft.fftshift(np.fft.fftn (np.fft.ifftshift(x),[Nx,Ny,Nz]))
ifft = lambda x: np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(x),[Nx,Ny,Nz]))

# functions for the projection 'G', and the product 'G : K^LT : (delta F)^T'
G      = lambda A2 : np.real( ifft( ddot42(Ghat4,fft(A2)) ) ).reshape(-1)
K_dF   = lambda dFm: trans2(ddot42(K4,trans2(dFm.reshape(3,3,Nx,Ny,Nz))))
G_K_dF = lambda dFm: G(K_dF(dFm))


# This is the CEIT constitutive model
import numpy as np



# material properties

E = 200.0e9
nu = 0.3
tau0 = 160.0e6
taus = 3031.3e6
C1 = 0.68
tausat = 140.0e6
tauv0 = 120.0e6
beta = 34.48
sigmaycem = 4200.0e6
b = 0.248
wtC = 0.77

props = [E,nu,tau0,taus,C1,tausat,tauv0,beta,sigmaycem,b,wtC]


# material history
colony_size = 10000
lamellar_spacing0 = 80.174169
lamellar_direction = [-0.21587,-0.938640,-0.268990] # direction cosines? are these the cos(angle) for each cardinal direction
in_plane_angle = 2.518623 # not sure what this is, but is said to be irrelevant in the document
effective_stress_ferrite = 0.0 # not required
eq_plastic_strain_ferrite = 0.0
hill_effective_stresses = [0,0] #not required
ferrite_stress = [0,0,0,0,0,0]
cem_stress = [0,0,0]
newton_it = 0 # not required
lamellar_direction0 = lamellar_direction
in_plane_angle0 = in_plane_angle
eq_plastic_strain_cem = 0 # not required

statOld = [
    colony_size, lamellar_spacing0, lamellar_direction[0], lamellar_direction[1],
    lamellar_direction[2], in_plane_angle, effective_stress_ferrite,
    eq_plastic_strain_ferrite, hill_effective_stresses[0],
    hill_effective_stresses[1],
    ferrite_stress[0],ferrite_stress[1],ferrite_stress[2],
    ferrite_stress[3],ferrite_stress[4],ferrite_stress[5],
    cem_stress[0],cem_stress[1],cem_stress[2],
    newton_it, lamellar_direction0[0],lamellar_direction0[1],
    lamellar_direction0[2], in_plane_angle0, eq_plastic_strain_cem
]


props = np.asarray(props)
statOld = np.asarray(statOld)

# storage for the state variables at every integration point
from random import gauss

e1 = np.array([1,0,0])
e2 = np.array([0,1,0])
e3 = np.array([0,0,1])


def make_rand_vector(dims):
    vec = [gauss(0, 1) for i in range(dims)]
    mag = sum(x**2 for x in vec) ** .5
    return [x/mag for x in vec]

stat_t = np.zeros([len(statOld),Nx,Ny,Nz])
for x in range(0,Nx):
    for y in range(0,Ny):
        for z in range(0,Nz):
            stat_t[:,x,y,z] = np.copy(statOld)

            V = make_rand_vector(3)
            magV = np.linalg.norm(V)
            
            stat_t[2:6,x,y,z] = [np.dot(V,e1)/magV, np.dot(V,e2)/magV, np.dot(V,e3)/magV, np.random.rand()*np.pi]
            stat_t[20:24,x,y,z] = stat_t[2:6,x,y,z]
            
                   


if inputPkl:
    print('Applying state variables from', inputPkl)
    stat_t = pickle.load(open(inputPkl,'rb'))
            
# have these as sorts of global variables
K4 = np.zeros([3,3,3,3,Nx,Ny,Nz])
P = np.zeros([3,3,Nx,Ny,Nz])

def stressCalculate(F_cell,Ft_cell,props,stat_t_cell):
    
    F_cell_reorder = F_cell
    Ft_cell_reorder = Ft_cell
    
    statNew,cauchyVoigt =         mesostress_v2.mesostressatintpoint_2_22(F_cell_reorder,Ft_cell_reorder,props,stat_t_cell)
    
    cauchyStress = np.zeros([3,3])
    cauchyStress[0,0] = cauchyVoigt[0]
    cauchyStress[1,1] = cauchyVoigt[1]
    cauchyStress[2,2] = cauchyVoigt[2]
    cauchyStress[0,1] = cauchyVoigt[3]
    cauchyStress[0,2] = cauchyVoigt[4]
    cauchyStress[1,2] = cauchyVoigt[5]
    cauchyStress[1,0] = cauchyStress[0,1]
    cauchyStress[2,0] = cauchyStress[0,2]
    cauchyStress[2,1] = cauchyStress[1,2]
    
    return statNew,cauchyStress


def constitutiveCEIT(F_local,F_t_local,stat_t_local):

    stat_new_local = np.copy(stat_t_local)
    P_local = np.zeros_like(F_local)
    cauchy_local = np.zeros_like(F_local)    
    stat_local = np.zeros_like(stat_t_local)    
    # This loop over every cell could be parallelised
    for x in range(0,Nx):
        for y in range(0,Ny):
            for z in range(0,Nz):
                
                F_cell = F_local[:,:,x,y,z]
                Ft_cell = F_t_local[:,:,x,y,z]
                stat_t_cell = np.copy(stat_t_local[:,x,y,z])

                statNew_cell,cauchy_stress = stressCalculate(F_cell,Ft_cell,props,stat_t_cell)

                J = np.linalg.norm(F_cell)
                
                Finv = np.linalg.inv(F_cell)
                
                P_cell = np.dot(Finv,J*cauchy_stress).transpose()
                
                P_local[:,:,x,y,z] = P_cell
                cauchy_local[:,:,x,y,z] = cauchy_stress
                
                stat_new_local[:,x,y,z] = np.copy(statNew_cell)
    
    return(P_local,stat_new_local,cauchy_local)


# In[12]:
def numericalTangent(F_local,F_t_local,stat_t_local):
    numK = np.zeros_like(K4)

    deltaF = 1e-2

    for k in range(0,3):
        for l in range(0,3):
            
            F1 = np.copy(F_local)
            F2 = np.copy(F_local)            

            for x in range(0,Nx):
                for y in range(0,Ny):
                    for z in range(0,Nz):            
                        F1[l,k,x,y,z] = F1[l,k,x,y,z] - deltaF
                        F2[l,k,x,y,z] = F2[l,k,x,y,z] + deltaF
            
            [P1,stat_new_1,cauchy_1] = constitutiveCEIT(F1,F_t_local,stat_t_local)
            [P2,stat_new_2,cauchy_2] = constitutiveCEIT(F2,F_t_local,stat_t_local)            
            for i in range(0,3):
                for j in range(0,3):
                    for x in range(0,Nx):
                        for y in range(0,Ny):
                            for z in range(0,Nz):            
                                numK[i,j,k,l,x,y,z] = (P2[j,i,x,y,z] - P1[j,i,x,y,z])/(2.0*deltaF)

    return numK

# phase indicator: square inclusion of volume fraction (3*3*15)/(11*13*15)
phase  = np.zeros([Nx,Ny,Nz]); phase[:3,:3,:] = 1.
# function to convert material parameters to grid of scalars
param  = lambda M0,M1: M0*np.ones([Nx,Ny,Nz])*(1.-phase)+                       M1*np.ones([Nx,Ny,Nz])*    phase

# material parameters, modified K and mu for pearlite
K      = param(166e9,166e9)  # bulk      modulus
mu     = param(77e9,77e9)  # shear     modulus
H      = param(0.004,0.008)  # hardening modulus
tauy0  = param(0.003,0.006)  # initial yield stress

# volume avg stress
def calculateAvgs(P_field,stat_field):
    cells = Nx*Ny*Nz
        
    P_avg = np.zeros([3,3])
    cauchy_avg = np.zeros([3,3])
    s_avg = 0
    
    for x in range(0,Nx):
        for y in range(0,Ny):
            for z in range(0,Nz):
                P_avg += P_field[:,:,x,y,z]/cells
                s_avg += stat_field[1,x,y,z]/cells
                cauchy_avg += cauchy[:,:,x,y,z]/cells
    return P_avg,s_avg,cauchy_avg

TIC = time.time()

print_freq = 1

# ---------------------------------- LOADING ----------------------------------


F_history = []
cauchy_stress_history = []
P_history = []
s_history = []
cauchy_history = []

# stress, deformation gradient, plastic strain, elastic Finger tensor
# NB "_t" signifies that it concerns the value at the previous increment
P      = np.zeros([3,3,Nx,Ny,Nz])
cauchy = np.zeros([3,3,Nx,Ny,Nz])
F      = np.array(I,copy=True)
F_t    = np.array(I,copy=True)

# initialize macroscopic incremental loading
ninc   = nIncrements
lam    = 0.0
barF   = np.array(I,copy=True)
barF_t = np.array(I,copy=True)

# the time increment
dt = T_final/ninc

# initial tangent operator: the elastic tangent
K4     = K*II+2.*mu*(I4s-1./3.*II)


t = 0


relaxationFactor = defaultRelaxationFactor


# create the output file 
outputFile = open(outputFile_loc,'w')
outputFile.write(
    "Time F_xx F_xy F_xz F_yx F_yy F_yz F_zx F_zy F_zz "
    "sigma_xx sigma_yy sigma_zz sigma_xy sigma_yz sigma_zx meanLamellarSpacing mean_eq_strain_ferr mean_eq_strain_cem\n"
    )

outputFile.write(
    '0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 ' +
    str(np.mean(stat_t[1,:,:,:])) + ' ' + str(np.mean(stat_t[7,:,:,:])) + ' ' +
    str(np.mean(stat_t[24,:,:,:]))
)



# incremental deformation
for inc in range(0,ninc):

    t += dt

    FNew = np.array([
        [F11(t),F12(t),F13(t)],
        [F21(t),F22(t),F23(t)],
        [F31(t),F32(t),F33(t)]
    ])
    
    
    #if (inc%print_freq == 0):
    print('=============================')
    print('inc: {0:d}'.format(inc))
        #try:
        #    print(P_history[-1])
        #except:
        #    pass

    barF      = np.array(I,copy=True)
    
    barF[0,0] = FNew[0,0]
    barF[0,1] = FNew[0,1]
    barF[0,2] = FNew[0,2]
    barF[1,0] = FNew[1,0]
    barF[1,1] = FNew[1,1]
    barF[1,2] = FNew[1,2]
    barF[2,0] = FNew[2,0]
    barF[2,1] = FNew[2,1]
    barF[2,2] = FNew[2,2]    
    
    #barF[0,0] =    (1.+lam)
    #barF[1,1] = 1./np.sqrt(barF[0,0])
    #barF[2,2] = barF[1,1]

    # store normalization
    Fn = np.linalg.norm(F)

    
    # first iteration residual: distribute "barF" over grid using "K4"
    b     = -G_K_dF(barF-barF_t)
    F    +=         barF-barF_t

    # parameters for Newton iterations: normalization and iteration counter
    Fn    = np.linalg.norm(F)
    iiter = 0

    # start with an approximation
    K4     = K*II+2.*mu*(I4s-1./3.*II)
    # try get away with the basic linear tangent for as
    # long as possible
    #usingNumericalGradient = False

    #relaxationFactor = defaultRelaxationFactor
    
    # iterate as long as the iterative update does not vanish   
    
    while True:

        print('$$$$$ iter:',iiter)
        if usingNumericalGradient:
            K4 = numericalTangent(F,F_t,stat_t)

        
        # solve linear system using the Conjugate Gradient iterative solver
        dFm,_ = sp.cg(tol=1.e-8,
          A   = sp.LinearOperator(shape=(F.size,F.size),matvec=G_K_dF,dtype='float'),
          b   = b,
        )

        # add solution of linear system to DOFs
        F += relaxationFactor*dFm.reshape(3,3,Nx,Ny,Nz)

        # compute residual stress and tangent, convert to residual
        P,stat_new,cauchy = constitutiveCEIT(F,F_t,stat_t)
        b          = -G(P)

        
        # check for convergence, print convergence info to screen
        residual = np.linalg.norm(dFm)/Fn
        if residual > 1:
            iiter = maxIts
        
        if (inc%print_freq == 0):
            print('{0:10.2e}'.format(residual))
        if np.linalg.norm(dFm)/Fn<1.e-5 and iiter>0: break

        # update Newton iteration counter
        iiter += 1   

        if iiter > maxIts:
            if relaxationFactor != defaultRelaxationFactor:
                raise Exception('Too many Newton iterations, exiting')
            else:
                print('Retrying increment with reduced relaxation')
                
                iiter = -1
                usingNumericalGradient = True
                F    = np.array(F_t   ,copy=True)
                # first iteration residual: distribute "barF" over grid using "K4"
                b     = -G_K_dF(barF-barF_t)
                F    +=         barF-barF_t

                # try a low relaxation value if it doesnt converge
                relaxationFactor = 0.1*defaultRelaxationFactor

    #print(barF[:20])
    TOC = time.time()

    time_taken = TOC - TIC
    print('Simulation time:',time_taken)
    
    # end-of-increment: update history
    barF_t = np.array(barF,copy=True)
    F_t    = np.array(F   ,copy=True)
    stat_t   = np.array(stat_new  ,copy=True)

    P_avg,s_avg,cauchy_avg = calculateAvgs(P,stat_new)
    #P_history.append(P_avg)
    #s_history.append(s_avg)
    #cauchy_history.append(cauchy_avg)

    outputFile.write(
        '\n' +
        str(t) + ' ' +
        str(FNew[0,0]) + ' ' + str(FNew[0,1]) + ' ' + str(FNew[0,2]) + ' ' +
        str(FNew[1,0]) + ' ' + str(FNew[1,1]) + ' ' + str(FNew[1,2]) + ' ' +
        str(FNew[2,0]) + ' ' + str(FNew[2,1]) + ' ' + str(FNew[2,2]) + ' ' +
        str(cauchy_avg[0,0]) + ' ' + str(cauchy_avg[1,1]) + ' ' +
        str(cauchy_avg[2,2]) + ' ' + str(cauchy_avg[0,1]) + ' ' +
        str(cauchy_avg[1,2]) + ' ' + str(cauchy_avg[2,1]) + ' ' +
        str(np.mean(stat_t[1,:,:,:])) + ' ' + str(np.mean(stat_t[7,:,:,:])) + ' ' +
        str(np.mean(stat_t[24,:,:,:]))
    )

outputFile.close()
# end 
