from coconut.tools import solver_available
from coconut.tests.solver_wrappers.abaqus_cse import abaqus_cse

import unittest

version = '2024'


@unittest.skipUnless(solver_available(f'abaqus_cse.v{version}'), f'abaqus_cse.v{version} not available')
class TestSolverWrapperAbaqusCSE2024Tube2D(abaqus_cse.TestSolverWrapperAbaqusCSETube2D):
    version = version
    setup_case = True

    def test_shear(self):
        super().test_shear()

    def test_repeat_iteration(self):
        super().test_repeat_iteration()

    def test_partitioning(self):
        super().test_partitioning()


@unittest.skipUnless(solver_available(f'abaqus_cse.v{version}'), f'abaqus_cse.v{version} not available')
class TestSolverWrapperAbaqusCSE2024Tube3D(abaqus_cse.TestSolverWrapperAbaqusCSETube3D):
    version = version
    setup_case = True

    def test_shear(self):
        super().test_shear()

    def test_repeat_iteration(self):
        super().test_repeat_iteration()

    def test_partitioning(self):
        super().test_partitioning()


if __name__ == '__main__':
    unittest.main()
