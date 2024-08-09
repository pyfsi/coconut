from coconut.tools import solver_available
from coconut.tests.solver_wrappers.openfoam import openfoam

import unittest

version_label = '10'


@unittest.skipUnless(solver_available(f'openfoam.v{version_label}'), f'openfoam.v{version_label} not available')
class TestSolverWrapperOpenFOAM10(openfoam.TestSolverWrapperOpenFOAM):
    version = version_label


if __name__ == '__main__':
    unittest.main()
