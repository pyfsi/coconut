from coconut.tools import solver_available
from coconut.tests.solver_wrappers.abaqus_line_load import abaqus

import unittest

version = '2024'


@unittest.skipUnless(solver_available(f'abaqus.line_load_v{version}'), f'abaqus.line_load_v{version} not available')
class TestSolverWrapperAbaqusLineLoad2024(abaqus.TestSolverWrapperAbaqusLineLoad):
    version = version
    setup_case = True


if __name__ == '__main__':
    unittest.main()
