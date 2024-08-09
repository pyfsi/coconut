from coconut.tools import create_instance, get_solver_env

import unittest
import numpy as np
import os
from os.path import join
import subprocess
import json
import shutil


class TestSolverWrapperAbaqusLineLoad(unittest.TestCase):
    version = None  # Abaqus version, e.g. 2022, set in subclass
    setup_case = True
    mp_name_in = 'YARN_load_points'
    mp_name_out = 'YARN_nodes'

    @classmethod
    def setUpClass(cls):
        dir_name = os.path.realpath(os.path.dirname(__file__))
        cls.file_name = join(dir_name, f'test_v{cls.version}/yarn/parameters.json')
        cls.working_dir = join(dir_name, f'test_v{cls.version}/yarn/CSM')
        dir_tmp = join(dir_name, f'test_v{cls.version}/yarn')

        # perform reference calculation
        with open(cls.file_name) as parameter_file:
            parameters = json.load(parameter_file)

        parameters['settings']['working_directory'] = os.path.relpath(cls.working_dir)  # set working directory
        solver_name = parameters['type'].replace('solver_wrappers.', '')
        env = get_solver_env(solver_name, dir_tmp)

        if cls.setup_case:
            p = subprocess.Popen('sh ' + join(dir_tmp, 'setup_abaqus.sh'), cwd=dir_tmp, shell=True, env=env)
            p.wait()

        # create the solver
        solver = create_instance(parameters)
        solver.initialize()
        interface_input = solver.get_interface_input()
        model_part = interface_input.get_model_part(cls.mp_name_in)

        # give value to variables
        traction = cls.get_line_load(model_part.x0)
        interface_input.set_variable_data(cls.mp_name_in, 'traction', traction)

        # step 1, coupling 1
        solver.initialize_solution_step()
        output1_1 = solver.solve_solution_step(interface_input)
        # step 1, coupling 2
        output1_2 = solver.solve_solution_step(interface_input)
        solver.finalize_solution_step()
        solver.output_solution_step()

        # save output for comparison, as input hasn't changed these should be the same
        cls.a1_1 = output1_1.get_variable_data(cls.mp_name_out, 'displacement')
        cls.a2_1 = output1_2.get_variable_data(cls.mp_name_out, 'displacement')

        # step 2 to 4
        for i in range(3):
            solver.initialize_solution_step()
            solver.solve_solution_step(interface_input)
            solver.finalize_solution_step()
            solver.output_solution_step()
        solver.finalize()

        # get data for solver without restart
        cls.interface_y_single_run = solver.get_interface_input()
        cls.interface_x_single_run = solver.get_interface_output()
        output_single_run = solver.get_interface_output()
        cls.a1 = output_single_run.get_variable_data(cls.mp_name_out, 'displacement')
        print(f'Max disp a1: {np.max(np.abs(cls.a1), axis=0)}')

    @classmethod
    def tearDownClass(cls):
        if cls.setup_case:
            shutil.rmtree(cls.working_dir)

    def setUp(self):
        with open(self.file_name) as parameter_file:
            self.parameters = json.load(parameter_file)
        self.parameters['settings']['working_directory'] = os.path.relpath(self.working_dir)  # set working directory

    @staticmethod
    def get_line_load(x):
        shear = np.zeros((x.shape[0], 3))
        shear[:, 1] = 1.5 * np.sin((x - x[-1]) / (x[0] - x[-1]) * np.pi)
        return shear

    def test_repeat_iteration(self):
        """
        Test whether repeating the same iteration yields the same results.
        """
        np.testing.assert_allclose(self.a2_1, self.a1_1, rtol=1e-15)

    def test_restart(self):
        """
        Test whether restarting at time step 2 and simulating 2 time steps yields the same displacement as when the
        simulation is run from time step 0 until time step 4.
        """

        # create solver which restarts at time step 2
        self.parameters['settings']['timestep_start'] = 2
        solver = create_instance(self.parameters)
        solver.initialize()
        interface_input = solver.get_interface_input()
        model_part = interface_input.get_model_part(self.mp_name_in)

        # give value to variables
        traction = self.get_line_load(model_part.x0)
        interface_input.set_variable_data(self.mp_name_in, 'traction', traction)

        # do step 3 and 4
        for i in range(2):
            solver.initialize_solution_step()
            solver.solve_solution_step(interface_input)
            solver.finalize_solution_step()
            solver.output_solution_step()
        solver.finalize()

        # compare output, as input hasn't changed these should be the same
        # get data for solver with restart
        interface_y_restart = solver.get_interface_input()
        interface_x_restart = solver.get_interface_output()
        output_restart = solver.get_interface_output()
        self.a3 = output_restart.get_variable_data(self.mp_name_out, 'displacement')
        print(f'\nMax disp a3: {np.max(np.abs(self.a3), axis=0)}')
        print(f'Max diff between a1 and a3: {np.abs(self.a1 - self.a3).max(axis=0)}')

        self.assertTrue(self.interface_y_single_run.has_same_model_parts(interface_y_restart))
        self.assertTrue(self.interface_x_single_run.has_same_model_parts(interface_x_restart))
        np.testing.assert_allclose(self.a3, self.a1, rtol=1e-10, atol=1e-17)

    def test_partitioning(self):
        """
        Test whether using 4 CPUs yields the same results as using a single one.
        """

        # adapt Parameters, create solver
        self.parameters['settings']['cores'] = 4
        solver = create_instance(self.parameters)
        solver.initialize()
        interface_input = solver.get_interface_input()
        model_part = interface_input.get_model_part(self.mp_name_in)

        # give value to variables
        traction = self.get_line_load(model_part.x0)
        interface_input.set_variable_data(self.mp_name_in, 'traction', traction)
        # do 4 steps
        for i in range(4):
            solver.initialize_solution_step()
            solver.solve_solution_step(interface_input)
            solver.finalize_solution_step()
            solver.output_solution_step()
        solver.finalize()

        # compare output, as input hasn't changed these should be the same
        output_4cores = solver.get_interface_output()
        self.a4 = output_4cores.get_variable_data(self.mp_name_out, 'displacement')
        print(f'\nMax disp a4: {np.max(np.abs(self.a4), axis=0)}')
        print(f'Max diff between a1 and a4: {np.abs(self.a1 - self.a4).max(axis=0)}')

        np.testing.assert_allclose(self.a4, self.a1, rtol=1e-10, atol=1e-17)


if __name__ == '__main__':
    unittest.main()
