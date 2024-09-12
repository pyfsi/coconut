from coconut.tools import solver_available
from coconut.tests.solver_wrappers.fluent_alm import fluent_alm

import unittest


version = '2024R2'


@unittest.skipUnless(solver_available(f'fluent_alm.v{version}'), f'fluent_alm.v{version} not available')
class TestSolverWrapperFluentALM2024R2Yarn(fluent_alm.TestSolverWrapperFluentALMYarn):
    version = version
    setup_case = True


if __name__ == '__main__':
    unittest.main()
