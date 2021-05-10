from coconut.tests.solver_wrappers.fluent import test_v2019R1 as TubeTest  # only from ... import ... works

import unittest


class TestSolverWrapperFluent2019R2Tube2D(TubeTest.TestSolverWrapperFluent2019R1Tube2D):
    version = '2019R2'
    setup_case = True


class TestSolverWrapperFluent2019R2Tube3D(TubeTest.TestSolverWrapperFluent2019R1Tube3D):
    version = '2019R2'
    setup_case = True


if __name__ == '__main__':
    unittest.main()
