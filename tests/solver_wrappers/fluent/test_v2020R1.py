from coconut.tests.solver_wrappers.fluent.fluent import TestSolverWrapperFluentTube2D, TestSolverWrapperFluentTube3D

import unittest


class TestSolverWrapperFluent2020R1Tube2D(TestSolverWrapperFluentTube2D):
    version = '2020R1'
    setup_case = True


class TestSolverWrapperFluent2020R1Tube3D(TestSolverWrapperFluentTube3D):
    version = '2020R1'
    setup_case = True


if __name__ == '__main__':
    unittest.main()
