from unittest import TestSuite
from coconut.tests.solver_wrappers.abaqus.Test_614 import TestSolverWrapperAbaqus614
from coconut.tests.solver_wrappers.fluent.test_v2019R1 import TestSolverWrapperFluent2019R1
from coconut.tests.solver_wrappers.fluent.test_v2019R2 import TestSolverWrapperFluent2019R2
from coconut.tests.solver_wrappers.fluent.test_v2019R3 import TestSolverWrapperFluent2019R3
from coconut.tests.solver_wrappers.fluent.test_v2020R1 import TestSolverWrapperFluent2020R1
from coconut.tests.solver_wrappers.openfoam.test_41 import TestSolverWrapperOpenFoam41
from coconut.tests.solver_wrappers.python.tube.test_tube_flow_solver import TestSolverWrapperTubeFlowSolver
from coconut.tests.solver_wrappers.python.tube.test_tube_ringmodel_solver import TestSolverWrapperTubeRingmodelSolver
from coconut.tests.solver_wrappers.python.tube.test_tube_structure_solver import TestSolverWrapperTubeStructureSolver

# test_cases = (TestSolverWrapperFluent2019R1, TestSolverWrapperFluent2019R2, TestSolverWrapperFluent2019R3, TestSolverWrapperFluent2020R1, TestSolverWrapperAbaqus614,
#               TestSolverWrapperOpenFoam41, TestSolverWrapperTubeFlowSolver, TestSolverWrapperTubeRingmodelSolver, TestSolverWrapperTubeStructureSolver)
test_cases = [TestSolverWrapperTubeStructureSolver,
              TestSolverWrapperTubeFlowSolver,
              TestSolverWrapperTubeRingmodelSolver]
               #TestSolverWrapperFluent2019R1]

def load_tests(loader, tests, pattern):
    suite = TestSuite()
    for test_class in test_cases:
        tests = loader.loadTestsFromTestCase(test_class)
        suite.addTests(tests)
    return suite