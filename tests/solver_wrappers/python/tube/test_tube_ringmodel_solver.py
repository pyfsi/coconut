from coconut.tools import create_instance

import unittest
import os
import json
import numpy as np
import subprocess


class TestSolverWrapperTubeRingmodelSolver(unittest.TestCase):

    def test_solver_wrapper_tube_ringmodel_solver(self):
        parameter_file_name = os.path.join(os.path.dirname(__file__), 'test_tube_ringmodel',
                                           'test_tube_ringmodel_solver.json')
        with open(parameter_file_name, 'r') as parameter_file:
            parameters = json.load(parameter_file)
        parameters_solver = parameters['solver_wrappers'][0]

        # if running from this folder
        if os.getcwd() == os.path.realpath(os.path.dirname(__file__)):
            parameters_solver['settings'].SetString('working_directory', 'test_tube_ringmodel/CSM')

        # setup case
        dir_tmp = os.path.join(os.path.realpath(os.path.dirname(__file__)), 'test_tube_ringmodel')
        p = subprocess.Popen(os.path.join(dir_tmp, 'setup_tube_ringmodel.sh'), cwd=dir_tmp, shell=True)
        p.wait()

        def get_dp(x):
            return 20 * np.sin(2 * np.pi / 0.05 * x)

        variable = 'pressure'
        model_part_name = 'wall'

        # test if same pressure always give same displacement
        if True:
            # create solver
            solver_1 = create_instance(parameters_solver)
            solver_2 = create_instance(parameters_solver)
            solvers = [solver_1, solver_2]
            for solver in solvers:
                solver.initialize()
                solver.initialize_solution_step()

            # change solver_1 to end pressure and solve
            interface_input = solver_1.get_interface_input()
            model_part = interface_input.get_model_part(model_part_name)
            data = [[get_dp(model_part.x0[i])] for i in range(model_part.size)]
            interface_input.set_variable_data(model_part_name, variable, np.array(data))
            output1_end = solver_1.solve_solution_step(interface_input).copy()

            # change solver_2 to intermediate pressure and solve
            interface_input = solver_2.get_interface_input()
            model_part = interface_input.get_model_part(model_part_name)
            data = [[0.5 * get_dp(model_part.x0[i])] for i in range(model_part.size)]
            interface_input.set_variable_data(model_part_name, variable, np.array(data))
            solver_2.solve_solution_step(interface_input).copy()

            # change solver_2 to end position and solve
            interface_input = solver_2.get_interface_input()
            model_part = interface_input.get_model_part(model_part_name)
            data = [[get_dp(model_part.x0[i])] for i in range(model_part.size)]
            interface_input.set_variable_data(model_part_name, variable, np.array(data))
            output2_end = solver_2.solve_solution_step(interface_input).copy()

            for solver in solvers:
                solver.finalize_solution_step()
                solver.finalize()

            # compare
            a1 = output1_end.get_interface_data()
            a2 = output2_end.get_interface_data()

            np.testing.assert_allclose(a1, a2, atol=1e-12)


if __name__ == '__main__':
    unittest.main()
