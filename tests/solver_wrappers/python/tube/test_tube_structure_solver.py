from coconut.tools import create_instance

import unittest
import os
import shutil
import json
import numpy as np


def get_dp(x):
    return 1500 * np.sin(2 * np.pi / 0.05 * x)


class TestSolverWrapperTubeStructureSolver(unittest.TestCase):

    def setUp(self):
        dir_name = os.path.realpath(os.path.dirname(__file__))  # path to python.tube directory

        # read settings
        parameter_file_name = os.path.join(dir_name, 'test_tube_structure/test_tube_structure_solver.json')
        with open(parameter_file_name, 'r') as parameter_file:
            self.parameters = json.load(parameter_file)
        self.variable = 'pressure'
        self.model_part_name = 'wall'

        # set working directory
        self.working_dir = os.path.join(dir_name, 'test_tube_structure/CSM')
        self.parameters['settings']['working_directory'] = os.path.relpath(self.working_dir)

        # setup
        shutil.rmtree(os.path.join(dir_name, self.working_dir), ignore_errors=True)
        os.mkdir(self.working_dir)
        shutil.copy(os.path.join(dir_name, 'test_tube_structure/setup_tube_structure/solver_parameters.json'),
                    self.working_dir)

    def test_pressure(self):
        # test if same pressure always give same displacement

        # create solver
        solver_1 = create_instance(self.parameters)
        solver_2 = create_instance(self.parameters)
        solvers = [solver_1, solver_2]
        for solver in solvers:
            solver.initialize()
            solver.initialize_solution_step()

        # change solver_1 to end pressure and solve
        interface_input = solver_1.get_interface_input()
        model_part = interface_input.get_model_part(self.model_part_name)
        data = [[get_dp(model_part.z0[i])] for i in range(model_part.size)]
        interface_input.set_variable_data(self.model_part_name, self.variable, np.array(data))
        output1_end = solver_1.solve_solution_step(interface_input).copy()

        # change solver_2 to intermediate pressure and solve
        interface_input = solver_2.get_interface_input()
        model_part = interface_input.get_model_part(self.model_part_name)
        data = [[0.5 * get_dp(model_part.z0[i])] for i in range(model_part.size)]
        interface_input.set_variable_data(self.model_part_name, self.variable, np.array(data))
        solver_2.solve_solution_step(interface_input).copy()

        # change solver_2 to end position and solve
        interface_input = solver_2.get_interface_input()
        model_part = interface_input.get_model_part(self.model_part_name)
        data = [[get_dp(model_part.z0[i])] for i in range(model_part.size)]
        interface_input.set_variable_data(self.model_part_name, self.variable, np.array(data))
        output2_end = solver_2.solve_solution_step(interface_input).copy()

        for solver in solvers:
            solver.finalize_solution_step()
            solver.output_solution_step()
            solver.finalize()

        # compare
        a1 = output1_end.get_interface_data()
        a2 = output2_end.get_interface_data()

        np.testing.assert_allclose(a1, a2, atol=1e-15)

    def tearDown(self):
        shutil.rmtree(self.working_dir)


if __name__ == '__main__':
    unittest.main()
