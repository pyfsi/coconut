import coconut.tests.solver_wrappers.fluent.test_v2019R1 as TubeTest

import unittest


class TestSolverWrapperFluent2019R3Tube2D(TubeTest.TestSolverWrapperFluent2019R1Tube2D):
    version = '2019R3'
    setup_case = True


class TestSolverWrapperFluent2019R3Tube3D(TubeTest.TestSolverWrapperFluent2019R1Tube3D):
    version = '2019R3'
    setup_case = True


if __name__ == '__main__':
    unittest.main()
