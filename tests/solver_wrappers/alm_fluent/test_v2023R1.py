from coconut.tools import solver_available
from coconut.tests.solver_wrappers.alm_fluent import alm_fluent

import unittest


version = '2023R1'


@unittest.skipUnless(solver_available(f'fluent.alm_v{version}'), f'fluent.alm_v{version} not available')
class TestSolverWrapperALMFluent2023R1Yarn(alm_fluent.TestSolverWrapperALMFluentYarn):
    version = version
    setup_case = True


if __name__ == '__main__':
    unittest.main()
