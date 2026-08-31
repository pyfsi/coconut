from coconut.tools import solver_available
from coconut.tests.solver_wrappers.openfoam_esi import openfoam_esi
import unittest

version_label = "2312"

is_solver_available = solver_available(f"openfoam_esi.v{version_label}")
error_msg = f'openfoam_esi.v{version_label} not available'
@unittest.skipUnless(is_solver_available, error_msg)
class TestSolverWrapperOpenFoamESI2312(openfoam_esi.TestSolverWrapperOpenFoamESI):
    version = version_label

if __name__=="__main__":
    unittest.main()