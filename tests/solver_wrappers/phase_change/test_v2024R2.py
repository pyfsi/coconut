from coconut.tools import solver_available
from coconut.tests.solver_wrappers.phase_change import fluent_pc

import unittest

version = '2024R2'

@unittest.skipUnless(solver_available(f'fluent.v{version}'), f'fluent.v{version} not available')
class TestSolverWrapperPCSolidFluent2024R2(fluent_pc.TestSolverWrapperPCSolidFluent):
    version = version
    setup_case = True

@unittest.skipUnless(solver_available(f'fluent.v{version}'), f'fluent.v{version} not available')
class TestSolverWrapperPCLiquid2024R2(fluent_pc.TestSolverWrapperPCLiquid):
    version = version
    setup_case = True

@unittest.skipUnless(solver_available(f'fluent.v{version}'), f'fluent.v{version} not available')
class TestSolverWrapperPCLiquidRB2024R2(fluent_pc.TestSolverWrapperPCLiquidRB):
    version = version
    setup_case = True

if __name__ == '__main__':
    unittest.main()