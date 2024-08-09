from coconut.tools import solver_available
from coconut.tests.solver_wrappers.abaqus import abaqus

import unittest

version = '2023'


@unittest.skipUnless(solver_available(f'abaqus.v{version}'), f'abaqus.v{version} not available')
class TestSolverWrapperAbaqus2023Tube2D(abaqus.TestSolverWrapperAbaqusTube2D):
    version = version
    setup_case = True


@unittest.skipUnless(solver_available(f'abaqus.v{version}'), f'abaqus.v{version} not available')
class TestSolverWrapperAbaqus2023Tube3D(abaqus.TestSolverWrapperAbaqusTube3D):
    version = version
    setup_case = True


if __name__ == '__main__':
    unittest.main()
