from coconut.tools import solver_available
from coconut.tests.solver_wrappers.kratos_structure import base_test_kratos_structure

import unittest

version = '60'


@unittest.skipUnless(solver_available(f'kratos_structure.v{version}'), f'kratos.structure_v{version} not available')
class TestSolverWrapperKratosStructure60(base_test_kratos_structure.BaseTestSolverWrapperKratosStructure):
    version_label = version


if __name__ == '__main__':
    unittest.main()
