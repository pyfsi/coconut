from coconut.tools import solver_available
from coconut.tests.solver_wrappers.openfoam_org import openfoam_org

import unittest

version_label = '10'


@unittest.skipUnless(solver_available(f'openfoam_org.v{version_label}'), f'openfoam_org.v{version_label} not available')
class TestSolverWrapperOpenFOAM10(openfoam_org.TestSolverWrapperOpenFoamORG):
    version = version_label


if __name__ == '__main__':
    unittest.main()
