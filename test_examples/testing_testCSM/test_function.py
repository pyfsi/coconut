#Import of utilities
import numpy as np

'''This is an example of test_function.py
To use test functions for the testing of a single solver, a file with this name should be included in the working directory.
The functions that need to be defined depend on the variables for the interface
The names of these functions are fixed as calculate_VARIABLE(x,y,z,n)
The functions receive the x, y and z-coordinate of the nodes in undeformed state and the current time step (n)
Several types of test can be grouped into this test_function.py file by creating additional classes
The name of the class to be used should be specified in the .json file containing the settings for the case
'''


class Simple_test():
    def calculate_DISPLACEMENT(self, x, y, z, n):
        """ Specify the displacement of a point on the interface relative to its original position"""
        if 0.01 < x < 0.024:
            disp = [0, 0, 0.001]
        else:
            disp = [0, 0, 0]
        return disp

    def calculate_PRESSURE(self, x, y, z, n):
        """ Specify the pressure on the surface"""
        if x > 0.01:
            pres = 1000
        else:
            pres = 0
        return pres

    def calculate_TRACTION(self, x, y, z, n):
        """ Specify the traction on the surface"""
        if x > 0.01:
            trac = [0, 0, 100]
        else:
            trac = [0, 0, 0]
        return trac


class Transient_test():
    def calculate_DISPLACEMENT(self, x, y, z, n):
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

    def calculate_PRESSURE(self, x, y, z, n):
        """ Specify the pressure on the surface"""
        if n < 5:
            if x > 0.01:
                pres = 1000
            else:
                pres = 0
        else:
            if x > 0.01:
                pres = -1000
            else:
                pres = 0
        return pres

    def calculate_TRACTION(self, x, y, z, n):
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
