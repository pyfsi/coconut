from unittest import TestSuite
# from coconut.tests.solver_wrappers.abaqus.Test_614 import TestSolverWrapperAbaqus614
from coconut.tests.solver_wrappers.fluent.test_v2019R1 import TestSolverWrapperFluent2019R1Tube2D, TestSolverWrapperFluent2019R1Tube3D
from coconut.tests.solver_wrappers.fluent.test_v2019R2 import TestSolverWrapperFluent2019R2Tube2D, TestSolverWrapperFluent2019R2Tube2D
from coconut.tests.solver_wrappers.fluent.test_v2019R3 import TestSolverWrapperFluent2019R3Tube2D, TestSolverWrapperFluent2019R3Tube2D
from coconut.tests.solver_wrappers.fluent.test_v2020R1 import TestSolverWrapperFluent2020R1Tube2D, TestSolverWrapperFluent2020R1Tube3D
# from coconut.tests.solver_wrappers.openfoam.test_41 import TestSolverWrapperOpenFoam41
from coconut.tests.solver_wrappers.python.tube.test_tube_flow_solver import TestSolverWrapperTubeFlowSolver
from coconut.tests.solver_wrappers.python.tube.test_tube_ringmodel_solver import TestSolverWrapperTubeRingmodelSolver
from coconut.tests.solver_wrappers.python.tube.test_tube_structure_solver import TestSolverWrapperTubeStructureSolver


test_cases = (
    # TestSolverWrapperAbaqus614,
    TestSolverWrapperFluent2019R1Tube2D, TestSolverWrapperFluent2019R1Tube3D,
    # TestSolverWrapperFluent2019R2Tube2D, TestSolverWrapperFluent2019R2Tube2D,
    # TestSolverWrapperFluent2019R3Tube2D, TestSolverWrapperFluent2019R3Tube2D,
    # TestSolverWrapperFluent2020R1Tube2D, TestSolverWrapperFluent2020R1Tube3D,
    # TestSolverWrapperOpenFoam41,
    TestSolverWrapperTubeFlowSolver,
    TestSolverWrapperTubeRingmodelSolver,
    TestSolverWrapperTubeStructureSolver
)


def load_tests(loader, tests, pattern):
    suite = TestSuite()
    for test_class in test_cases:
        tests = loader.loadTestsFromTestCase(test_class)
        suite.addTests(tests)
    return suite