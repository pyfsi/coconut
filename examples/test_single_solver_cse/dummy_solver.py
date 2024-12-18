import pickle

import matplotlib.pyplot as plt
import numpy as np
from coconut import tools
from scipy import interpolate

"""
This is an example of dummy_solver.py.
To use test functions for the testing of a single solver,
a file with this name should be included in the working directory.
There are two ways to create valid test classes:

1. With a method get_input(interface_input_to, n), setting this interface to the desired value.
To this end, its useful to use the optional initialize(interface_input_to, solver_index) method,
which is called at the start of the simulation and where interface_input_to is the input interface of tested solver.

2. With methods named calculate_<variable>(x, y, z, n), where <variable> is the variable on the interface.
The functions receive the x, y and z-coordinate of the data locations in undeformed state and the current time step (n).
They have to return a list or numpy array of 1 or 3 elements for a scalar or vector, respectively.

Several types of test can be grouped into this dummy_solver.py file by creating additional classes.
The name of the class to be used should be specified in the .json file containing the settings for the case.
"""


class TransientTest:
    def calculate_pressure(self, x, y, z, n):
        """ Specify the pressure on the surface"""
        return [10]

    def calculate_traction(self, x, y, z, n):
        """ Specify the traction on the surface"""
        return [1, 2, 3]
