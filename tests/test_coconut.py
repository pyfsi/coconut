from coconut.data_structure import KratosUnittest

from coconut.tests.coupled_solvers.test_gauss_seidel import TestCoupledSolverGaussSeidel
from coconut.tests.coupled_solvers.test_relaxation import TestCoupledSolverRelaxation
from coconut.tests.coupled_solvers.test_aitken import TestCoupledSolverAitken
from coconut.tests.coupled_solvers.test_iqni import TestCoupledSolverIQNI
from coconut.tests.coupled_solvers.test_ibqn import TestCoupledSolverIBQN
from coconut.tests.coupled_solvers.models.test_ls import TestModelLS
from coconut.tests.coupled_solvers.models.test_mv import TestModelMV

from coconut.tests.convergence_criteria.test_absolute_norm import TestConvergenceCriterionAbsoluteNorm
from coconut.tests.convergence_criteria.test_and import TestConvergenceCriterionAnd
from coconut.tests.convergence_criteria.test_iteration_limit import TestConvergenceCriterionIterationLimit
from coconut.tests.convergence_criteria.test_or import TestConvergenceCriterionOr
from coconut.tests.convergence_criteria.test_relative_norm import TestConvergenceCriterionRelativeNorm

from coconut.tests.mappers.test_interpolator import TestMapperInterpolator
from coconut.tests.mappers.test_nearest import TestMapperNearest
from coconut.tests.mappers.test_linear import TestMapperLinear
from coconut.tests.mappers.test_radial_basis import TestMapperRadialBasis
from coconut.tests.mappers.test_permutation import TestMapperPermutation
from coconut.tests.mappers.test_axisymmetric_2d_to_3d import TestMapperAxisymmetric2DTo3D
from coconut.tests.mappers.test_axisymmetric_3d_to_2d import TestMapperAxisymmetric3DTo2D
from coconut.tests.mappers.test_combined import TestMapperCombined

from coconut.tests.predictors.test_predictor import TestPredictor
from coconut.tests.predictors.test_constant import TestPredictorConstant
from coconut.tests.predictors.test_linear import TestPredictorLinear
from coconut.tests.predictors.test_legacy import TestPredictorLegacy
from coconut.tests.predictors.test_quadratic import TestPredictorQuadratic
from coconut.tests.predictors.test_cubic import TestPredictorCubic

from coconut.tests.solver_wrappers.python.tube.test_tube_flow_solver import TestSolverWrapperTubeFlowSolver
from coconut.tests.solver_wrappers.python.tube.test_tube_ringmodel_solver import TestSolverWrapperTubeRingmodelSolver
from coconut.tests.solver_wrappers.python.tube.test_tube_structure_solver import TestSolverWrapperTubeStructureSolver
from coconut.tests.solver_wrappers.fluent.test_v2019R1 import TestSolverWrapperFluent2019R1
from coconut.tests.solver_wrappers.fluent.test_v2019R2 import TestSolverWrapperFluent2019R2
from coconut.tests.solver_wrappers.fluent.test_v2019R3 import TestSolverWrapperFluent2019R3
from coconut.tests.solver_wrappers.fluent.test_v2020R1 import TestSolverWrapperFluent2020R1
from coconut.tests.solver_wrappers.abaqus.test_v614 import TestSolverWrapperAbaqus614
from coconut.tests.solver_wrappers.openfoam.test_41 import TestSolverWrapperOpenFoam41

from coconut.tests.data_structure.test_parameters import TestPyKratosParameters
from coconut.tests.data_structure.test_variables import TestPyKratosVariables
from coconut.tests.data_structure.test_cosimulation_interface import TestInterface


def AssembleTestSuites():
    ''' Populates the test suites to run.

    Populates the test suites to run. At least, it should populate the suites:
    "small", "nightly" and "all"

    Return
    ------

    suites: A dictionary of suites
        The set of suites with its test_cases added.
    '''

    suites = KratosUnittest.KratosSuites

    smallSuite = suites['small']  # These tests are executed by the continuous integration tool

    smallSuite.addTest(TestCoupledSolverGaussSeidel("test_coupled_solver_gauss_seidel"))
    smallSuite.addTest(TestCoupledSolverRelaxation("test_coupled_solver_relaxation"))
    smallSuite.addTest(TestCoupledSolverAitken("test_coupled_solver_aitken"))
    smallSuite.addTest(TestCoupledSolverIQNI("test_coupled_solver_iqni"))
    smallSuite.addTest(TestCoupledSolverIBQN("test_coupled_solver_ibqn"))
    smallSuite.addTest(TestModelLS("test_model_ls"))
    smallSuite.addTest(TestModelMV("test_model_mv"))

    smallSuite.addTest(TestConvergenceCriterionAbsoluteNorm("test_convergence_criterion_absolute_norm"))
    smallSuite.addTest(TestConvergenceCriterionAnd("test_convergence_criterion_and"))
    smallSuite.addTest(TestConvergenceCriterionIterationLimit("test_convergence_criterion_iteration_limit"))
    smallSuite.addTest(TestConvergenceCriterionOr("test_convergence_criterion_or"))
    smallSuite.addTest(TestConvergenceCriterionRelativeNorm("test_convergence_criterion_relative_norm"))

    smallSuite.addTest(TestMapperInterpolator("test_mapper_interpolator"))
    smallSuite.addTest(TestMapperNearest("test_mapper_nearest"))
    smallSuite.addTest(TestMapperLinear("test_mapper_linear"))
    smallSuite.addTest(TestMapperRadialBasis("test_mapper_radial_basis"))
    smallSuite.addTest(TestMapperPermutation("test_mapper_permutation"))
    smallSuite.addTest(TestMapperAxisymmetric2DTo3D("test_mapper_axisymmetric_2d_to_3d"))
    smallSuite.addTest(TestMapperAxisymmetric3DTo2D("test_mapper_axisymmetric_3d_to_2d"))
    smallSuite.addTest(TestMapperCombined("test_mapper_combined"))

    smallSuite.addTest(TestPredictor("test_predictor"))
    smallSuite.addTest(TestPredictorConstant("test_predictor_constant"))
    smallSuite.addTest(TestPredictorLinear("test_predictor_linear"))
    smallSuite.addTest(TestPredictorLegacy("test_predictor_legacy"))
    smallSuite.addTest(TestPredictorQuadratic("test_predictor_quadratic"))
    smallSuite.addTest(TestPredictorCubic("test_predictor_cubic"))

    smallSuite.addTest(TestSolverWrapperTubeFlowSolver("test_solver_wrapper_tube_flow_solver"))
    smallSuite.addTest(TestSolverWrapperTubeRingmodelSolver("test_solver_wrapper_tube_ringmodel_solver"))
    smallSuite.addTest(TestSolverWrapperTubeStructureSolver("test_solver_wrapper_tube_structure_solver"))
    # smallSuite.addTest(TestSolverWrapperFluent2019R1("test_solver_wrapper_fluent_2019R1"))  # duration ~500s
    # smallSuite.addTest(TestSolverWrapperFluent2019R2("test_solver_wrapper_fluent_2019R2"))
    # smallSuite.addTest(TestSolverWrapperFluent2019R3("test_solver_wrapper_fluent_2019R3"))
    # smallSuite.addTest(TestSolverWrapperFluent2020R1("test_solver_wrapper_fluent_2020R1"))
    # smallSuite.addTest(TestSolverWrapperAbaqus614("test_solver_wrapper_abaqus_614"))  # duration ~500s
    # smallSuite.addTest(TestSolverWrapperOpenFoam41('test_model_part_nodes_with_different_cores'))
    # smallSuite.addTest(TestSolverWrapperOpenFoam41('test_displacement_on_nodes'))
    # smallSuite.addTest(TestSolverWrapperOpenFoam41('test_displacement_on_nodes'))
    #
    smallSuite.addTest(TestPyKratosParameters("test_pykratos_parameters"))
    smallSuite.addTest(TestPyKratosVariables("test_pykratos_variables"))
    smallSuite.addTest(TestInterface("test_cosimulation_interface"))

    nightlySuite = suites['nightly']  # These tests are executed in the nightly build
    nightlySuite.addTests(smallSuite)

    validationSuite = suites['validation']   # These tests are very long and should not be in nightly, for validation

    # Create a test suit that contains all tests:
    allSuite = suites['all']
    allSuite.addTests(nightlySuite)  # Already contains the smallSuite
    allSuite.addTests(validationSuite)

    return suites


if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())
