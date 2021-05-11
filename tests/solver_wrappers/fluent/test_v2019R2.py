from coconut.tests.solver_wrappers.fluent.fluent import TestSolverWrapperFluentTube2D, TestSolverWrapperFluentTube3D

import unittest


class TestSolverWrapperFluent2019R2Tube2D(TestSolverWrapperFluentTube2D):
    version = '2019R2'
    setup_case = True


class TestSolverWrapperFluent2019R2Tube3D(TestSolverWrapperFluentTube3D):
    version = '2019R2'
    setup_case = True


if __name__ == '__main__':
    unittest.main()
