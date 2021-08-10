from coconut.tools import solver_available
from coconut.tests.solver_wrappers.openfoam_extend import openfoam_extend

import unittest

version_label = '41'


@unittest.skipUnless(solver_available(f'openfoam_extend.v{version_label}'), f'openfoam_extend.v{version_label} not available')
class TestSolverWrapperOpenFOAMExtend41(openfoam_extend.TestSolverWrapperOpenFOAMExtend):
    version = version_label


if __name__ == "__main__":
    unittest.main()
