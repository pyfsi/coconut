from coconut.tools import solver_available
from coconut.tests.solver_wrappers.fluent_alm import fluent_alm

import unittest


version = '2025R2'


@unittest.skipUnless(solver_available(f'fluent_alm.v{version}'), f'fluent_alm.v{version} not available')
class TestSolverWrapperFluentALM2025R2Yarn(fluent_alm.TestSolverWrapperFluentALMYarn):
    version = version
    setup_case = True


if __name__ == '__main__':
    unittest.main()
