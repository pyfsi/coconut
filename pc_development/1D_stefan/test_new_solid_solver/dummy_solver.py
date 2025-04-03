import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt

"""
This is an example of dummy_solver.py for thermal development.
To use test functions for the testing of a single solver, a file with this name should be included in the working directory.
The functions that need to be defined depend on the variables for the interface.
The names of these functions are fixed as calculate_<variable>(x,y,z,n).
The functions receive the x, y and z-coordinate of the nodes in undeformed state and the current time step (n).
They have to return a list or numpy array of 1 or 3 elements for a scalar or vector, respectively.
Several types of test can be grouped into this dummy_solver.py file by creating additional classes.
The name of the class to be used should be specified in the .json file containing the settings for the case.
"""


class SolidTest:
    def calculate_temperature(self, x, y, z, n):
        temp = [299.15]
        return temp

    def calculate_heat_flux(self, x, y, z, n):
        #hf = [-300 if n == 1 else 0]
        hf = [-300]
        return hf

    def calculate_displacement(self, x, y, z, n):
        dt = 1.0 # s
        Q = 300 # W/m^2
        rho = 870 # kg/m^3
        LH = 179000 # J/kg
        dx = -Q*dt/(rho*LH) # m
        disp = [0 if n == 1 else dx, 0, 0]
        return disp

class FluidTest:
    def calculate_temperature(self, x, y, z, n):
        temp = [299.15]
        return temp

    def calculate_heat_flux(self, x, y, z, n):
        hf = [0]
        return hf

    def calculate_displacement(self, x, y, z, n):
        disp = [-0.001*n, 0, 0]
        return disp