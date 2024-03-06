from coconut.tools import solver_available
from coconut.tests.solver_wrappers.kratos_structure import kratos_structure

import unittest

version = '94'


@unittest.skipUnless(solver_available(f'kratos_structure.v{version}'), f'kratos.structure_v{version} not available')
class TestSolverWrapperKratosStructure94(kratos_structure.BaseTestSolverWrapperKratosStructure):
    version_label = version


if __name__ == '__main__':
    unittest.main()
