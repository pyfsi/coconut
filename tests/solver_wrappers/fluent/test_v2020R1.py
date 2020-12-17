from coconut.tests.solver_wrappers.fluent.test_v2019R1 \
    import TestSolverWrapperFluent2019R1Tube2D, TestSolverWrapperFluent2019R1Tube3D

import unittest


class TestSolverWrapperFluent2020R1Tube2D(TestSolverWrapperFluent2019R1Tube2D):
    version = '2020R1'
    setup_case = True


class TestSolverWrapperFluent2020R1Tube3D(TestSolverWrapperFluent2019R1Tube3D):
    version = '2020R1'
    setup_case = True


if __name__ == '__main__':
    unittest.main()
