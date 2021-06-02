from coconut.tools import solver_available
from coconut.tests.solver_wrappers.openfoam import openfoam

import unittest


version = '8'


@unittest.skipUnless(solver_available(f'openfoam.v{version}'), f'openfoam.v{version} not available')
class TestSolverWrapperOpenFOAM8(openfoam.TestSolverWrapperOpenFOAM):
    version = version


if __name__ == '__main__':
    unittest.main()
