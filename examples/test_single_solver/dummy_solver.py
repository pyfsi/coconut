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
    def calculate_displacement(self, x, y, z, n):
        """ Specify the displacement of a point on the interface relative to its original position"""
        if n < 5:
            if 0.01 < x < 0.024:
                disp = [0, 0, 0.001]
            else:
                disp = [0, 0, 0]
        else:
            if 0.01 < x < 0.024:
                disp = [0, 0, -0.001]
            else:
                disp = [0, 0, 0]
        return disp

    def calculate_pressure(self, x, y, z, n):
        """ Specify the pressure on the surface"""
        if n < 5:
            if x > 0.01:
                pres = [1000]
            else:
                pres = [0]
        else:
            if x > 0.01:
                pres = [-1000]
            else:
                pres = [0]
        return pres

    def calculate_traction(self, x, y, z, n):
        """ Specify the traction on the surface"""
        if n < 5:
            if x > 0.01:
                trac = [0, 0, 100]
            else:
                trac = [0, 0, 0]
        else:
            if x > 0.01:
                trac = [0, 100, 0]
            else:
                trac = [0, 0, 0]
        return trac


class PythonSolverTest:
    def calculate_displacement(self, x, y, z, n):
        """ Specify the displacement of a point on the interface relative to its original position"""
        if n < 5:
            if 0.01 < z < 0.024:
                disp = [0, 0.0002, 0]
            else:
                disp = [0, 0, 0]
        else:
            if 0.01 < z < 0.024:
                disp = [0, -0.0001, 0]
            else:
                disp = [0, 0, 0]
        return disp

    def calculate_pressure(self, x, y, z, n):
        """ Specify the pressure on the surface"""
        if n < 5:
            if z > 0.01:
                pres = [1000]
            else:
                pres = [0]
        else:
            if z > 0.01:
                pres = [-1000]
            else:
                pres = [0]
        return pres

    def calculate_traction(self, x, y, z, n):
        """ Specify the traction on the surface"""
        if n < 5:
            if z > 0.01:
                trac = [0, 0, 100]
            else:
                trac = [0, 0, 0]
        else:
            if z > 0.01:
                trac = [0, 100, 0]
            else:
                trac = [0, 0, 0]
        return trac


# assumes time step size doesn't change!
class PickleData:
    # path to the pickle file used for input starting from the same location as this file
    input_pickle_file_name = 'input_case_results.pickle'
    # parameters for the mapper used to map the input pickle data to the interface of the solver
    mapper_parameters = {
        "type": "mappers.interface",
        "settings": {
            "type": "mappers.radial_basis",
            "settings": {
                "directions": [
                    "x",
                    "y",
                    "z"
                ]
            }
        }
    }

    def __init__(self):
        self.mapper = tools.create_instance(self.mapper_parameters)
        self.interface_input_from = None
        self.solution = None

    def initialize(self, interface_input_to, solver_index):
        with open(self.input_pickle_file_name, 'rb') as f:
            data = pickle.load(f)
        xy = 'x' if solver_index == 0 else 'y'
        self.interface_input_from = data[f'interface_{xy}']
        self.solution = data[f'solution_{xy}']
        self.mapper.initialize(self.interface_input_from, interface_input_to)

    def get_input(self, interface_input_to, n):
        self.interface_input_from.set_interface_data(self.solution[:, n])
        self.mapper(self.interface_input_from, interface_input_to)


class InterpolatedData:
    def __init__(self):
        abscissa = np.array([-0.025, -0.02433333, -0.02366667, -0.023, -0.02233333,
                             -0.02166667, -0.021, -0.02, -0.01666667, -0.01333333,
                             -0.01, -0.00666667, -0.00333333, 0., 0.00333333,
                             0.00666667, 0.01, 0.01333333, 0.01666667, 0.02,
                             0.02066667, 0.02133333, 0.022, 0.02266667, 0.02333333,
                             0.024, 0.02466667])
        y_displacement = np.array([1.16621303e-10, 1.95088490e-04, 3.60966621e-04, 4.83640889e-04,
                                   5.66206281e-04, 6.17896253e-04, 6.48283040e-04, 6.70057566e-04,
                                   6.78144584e-04, 6.75857193e-04, 6.74728085e-04, 6.73872500e-04,
                                   6.73091675e-04, 6.72355520e-04, 6.71655600e-04, 6.70986983e-04,
                                   6.70345131e-04, 6.69822764e-04, 6.70441553e-04, 6.63535153e-04,
                                   6.52518454e-04, 6.31071264e-04, 5.92546400e-04, 5.28036801e-04,
                                   4.27512507e-04, 2.83584753e-04, 9.96992219e-05])
        self.interpolant = interpolate.splrep(abscissa, y_displacement, s=0)
        plt.plot(abscissa, y_displacement)

    def calculate_displacement(self, x, y, z, n):
        """ Specify the displacement of a point on the interface relative to its original position"""
        if -0.025 <= x <= 0.025:
            disp = [0, float(interpolate.splev(x, self.interpolant, der=0)), 0]
        else:
            disp = [0, 0, 0]
        return disp
