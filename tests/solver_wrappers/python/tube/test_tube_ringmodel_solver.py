from coconut import data_structure
from coconut.data_structure import KratosUnittest
from coconut.coupling_components.tools import CreateInstance

import numpy as np
import os
import subprocess


class TestSolverWrapperTubeRingmodelSolver(KratosUnittest.TestCase):
    def assertArrayAlmostEqual(self, a1, a2, delta=None):
        ls1 = list(a1)
        ls2 = list(a2)
        try:
            self.assertEqual(ls1, ls2)
        except AssertionError:
            for i in range(len(ls1)):
                self.assertAlmostEqual(ls1[i], ls2[i], delta=delta)

    def test_solver_wrapper_tube_ringmodel_solver(self):
        parameter_file_name = os.path.join(os.path.dirname(__file__), 'test_tube_ringmodel',
                                           'test_tube_ringmodel_solver.json')
        with open(parameter_file_name, 'r') as parameter_file:
            parameters = data_structure.Parameters(parameter_file.read())
        parameters_solver = parameters['solver_wrappers'][0]

        # if running from this folder
        if os.getcwd() == os.path.realpath(os.path.dirname(__file__)):
            parameters_solver['settings'].SetString('working_directory', 'test_tube_ringmodel/CSM')

        # "global" definitions
        pressure = vars(data_structure)['PRESSURE']

        # setup case
        dir_tmp = os.path.join(os.path.realpath(os.path.dirname(__file__)), 'test_tube_ringmodel')
        p = subprocess.Popen(os.path.join(dir_tmp, 'setup_tube_ringmodel.sh'), cwd=dir_tmp, shell=True)
        p.wait()

        def get_dp(x):
            return 20 * np.sin(2 * np.pi / 0.05 * x)

        # test if same pressure always give same displacement
        if True:
            # create solver
            solver_1 = CreateInstance(parameters_solver)
            solver_2 = CreateInstance(parameters_solver)
            solvers = [solver_1, solver_2]
            for solver in solvers:
                solver.Initialize()
                solver.InitializeSolutionStep()

            # change solver_1 to end pressure and solve
            mp = solver_1.model['wall']
            for node in mp.Nodes:
                node.SetSolutionStepValue(pressure, 0, get_dp(node.X0))
            output1_end = solver_1.SolveSolutionStep(solver_1.GetInterfaceInput()).deepcopy()

            # change solver_2 to intermediate pressure and solve
            for node in mp.Nodes:
                node.SetSolutionStepValue(pressure, 0, 0.5 * get_dp(node.X0))
            solver_2.SolveSolutionStep(solver_2.GetInterfaceInput()).deepcopy()

            # change solver_2 to end pressure and solve
            for node in mp.Nodes:
                node.SetSolutionStepValue(pressure, 0, get_dp(node.X0))
            output2_end = solver_2.SolveSolutionStep(solver_2.GetInterfaceInput()).deepcopy()

            for solver in solvers:
                solver.FinalizeSolutionStep()
                solver.Finalize()

            # compare
            a1 = output1_end.GetNumpyArray()
            a2 = output2_end.GetNumpyArray()

            self.assertArrayAlmostEqual(a1, a2, delta=1e-12)


if __name__ == '__main__':
    KratosUnittest.main()
