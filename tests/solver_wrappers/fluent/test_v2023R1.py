from coconut.tools import solver_available
from coconut.tests.solver_wrappers.fluent import fluent

import unittest


version = '2023R1'


@unittest.skipUnless(solver_available(f'fluent.v{version}'), f'fluent.v{version} not available')
class TestSolverWrapperFluent2023R1Tube2D(fluent.TestSolverWrapperFluentTube2D):
    version = version
    setup_case = True


@unittest.skipUnless(solver_available(f'fluent.v{version}'), f'fluent.v{version} not available')
class TestSolverWrapperFluent2023R1Tube3D(fluent.TestSolverWrapperFluentTube3D):
    version = version
    setup_case = True


if __name__ == '__main__':
    unittest.main()
