from coconut.tests.solver_wrappers.fluent import fluent

import unittest


class TestSolverWrapperFluent2019R3Tube2D(fluent.TestSolverWrapperFluentTube2D):
    version = '2019R3'
    setup_case = True


class TestSolverWrapperFluent2019R3Tube3D(fluent.TestSolverWrapperFluentTube3D):
    version = '2019R3'
    setup_case = True


if __name__ == '__main__':
    unittest.main()
