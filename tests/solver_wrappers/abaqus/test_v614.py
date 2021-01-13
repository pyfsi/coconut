from coconut.coupling_components.tools import create_instance

import unittest
import numpy as np
import os
from os.path import join
import subprocess
import json


class TestSolverWrapperAbaqus614Tube2D(unittest.TestCase):
    setup_case = True
    dir_tmp = join(os.path.realpath(os.path.dirname(__file__)), 'test_v614/tube2d')
    file_name = join(os.path.dirname(__file__), 'test_v614/tube2d/parameters.json')
    p = 1500
    shear = np.array([0, 0, 0])
    mp_name_in = 'BEAMINSIDEMOVING_load_points'
    mp_name_out = 'BEAMINSIDEMOVING_nodes'

    @classmethod
    def setUpClass(cls):
        if cls.setup_case:
            p_setup_abaqus = subprocess.Popen(os.path.join(cls.dir_tmp, 'setup_abaqus.sh'), cwd=cls.dir_tmp, shell=True)
            p_setup_abaqus.wait()

        # perform reference calculation
        with open(cls.file_name, 'r') as parameter_file:
            parameters = json.load(parameter_file)

        print(f"Working directory reference: {parameters['settings']['working_directory']}")
        # create the solver
        solver = create_instance(parameters)
        interface_input = solver.get_interface_input()

        # give value to variables
        pressure = interface_input.get_variable_data(cls.mp_name_in, 'pressure')
        pressure[:] = cls.p
        interface_input.set_variable_data(cls.mp_name_in, 'pressure', pressure)
        traction = interface_input.get_variable_data(cls.mp_name_in, 'traction')
        traction[:, :] = cls.shear
        interface_input.set_variable_data(cls.mp_name_in, 'traction', traction)

        solver.initialize()
        # step 1, coupling 1
        solver.initialize_solution_step()
        output1_1 = solver.solve_solution_step(interface_input)
        # step 1, coupling 2
        output1_2 = solver.solve_solution_step(interface_input)
        solver.finalize_solution_step()

        # compare output, as input hasn't changed these should be the same
        a1 = output1_1.get_variable_data(cls.mp_name_out, 'displacement').copy()
        a2 = output1_2.get_variable_data(cls.mp_name_out, 'displacement').copy()

        # normalize data and compare
        mean = np.mean(a1, axis=0)
        ref = np.abs(a1 - mean).max()
        a1n = (a1 - mean) / ref
        a2n = (a2 - mean) / ref
        np.testing.assert_allclose(a1n, a2n, atol=1e-12)

        # step 2 to 4
        for i in range(3):
            solver.initialize_solution_step()
            solver.solve_solution_step(interface_input)
            solver.finalize_solution_step()
        solver.finalize()

        # get data for solver without restart
        output_single_run = solver.get_interface_output()
        cls.a1 = output_single_run.get_variable_data(cls.mp_name_out, 'displacement').copy()
        cls.mean = np.mean(cls.a1, axis=0)
        cls.ref = np.abs(cls.a1 - cls.mean).max()
        cls.a1n = (cls.a1 - cls.mean) / cls.ref
        cls.mean_disp_y_no_shear = np.mean(cls.a1[:, 1])/cls.a1.shape[0]
        cls.mean_disp_x_no_shear = np.mean(cls.a1[:, 0])/cls.a1.shape[0]   # needed for 3D test case

    def setUp(self):
        with open(self.file_name, 'r') as parameter_file:
            self.parameters = json.load(parameter_file)

    def test_restart(self):
        # test if restart option works correctly
        print(f"Working directory for test_restart: {self.parameters['settings']['working_directory']}")

        # create solver which restarts at time step 2
        self.parameters['settings']['timestep_start'] = 2
        solver = create_instance(self.parameters)
        interface_input = solver.get_interface_input()

        # give value to variables
        pressure = interface_input.get_variable_data(self.mp_name_in, 'pressure')
        pressure[:] = self.p
        interface_input.set_variable_data(self.mp_name_in, 'pressure', pressure)
        traction = interface_input.get_variable_data(self.mp_name_in, 'traction')
        traction[:, :] = self.shear
        interface_input.set_variable_data(self.mp_name_in, 'traction', traction)

        solver.initialize()
        # do step 3 and 4
        for i in range(2):
            solver.initialize_solution_step()
            solver.solve_solution_step(interface_input)
            solver.finalize_solution_step()
        solver.finalize()

        # compare output, as input hasn't changed these should be the same
        # get data for solver with restart
        output_restart = solver.get_interface_output()
        self.a3 = output_restart.get_variable_data(self.mp_name_out, 'displacement').copy()

        # normalize data and compare
        self.a3n = (self.a3 - self.mean) / self.ref
        print(f"Max diff between a1n and a3n: {np.abs(self.a1n - self.a3n).max()} \n")

        np.testing.assert_allclose(self.a1n, self.a3n, atol=1e-12)  # correct to consider absolute tolerance only?

    def test_partitioning(self):
        # test whether using 4 CPUs gives the same results as using a single one
        print(f"Working directory for test_partitioning: {self.parameters['settings']['working_directory']}")

        # adapt Parameters, create solver
        self.parameters['settings']['cores'] = 4
        solver = create_instance(self.parameters)
        interface_input = solver.get_interface_input()

        # give value to variables
        pressure = interface_input.get_variable_data(self.mp_name_in, 'pressure')
        pressure[:] = self.p
        interface_input.set_variable_data(self.mp_name_in, 'pressure', pressure)
        traction = interface_input.get_variable_data(self.mp_name_in, 'traction')
        traction[:, :] = self.shear
        interface_input.set_variable_data(self.mp_name_in, 'traction', traction)

        # do 4 steps
        solver.initialize()
        for i in range(4):
            solver.initialize_solution_step()
            solver.solve_solution_step(interface_input)
            solver.finalize_solution_step()
        solver.finalize()

        # compare output, as input hasn't changed these should be the same
        output_4cores = solver.get_interface_output()
        self.a4 = output_4cores.get_variable_data(self.mp_name_out, 'displacement').copy()

        # normalize data and compare
        self.a4n = (self.a4 - self.mean) / self.ref
        print(f"Max diff between a1n and a4n: {np.abs(self.a1n - self.a4n).max()} \n")

        np.testing.assert_allclose(self.a4n, self.a1n, atol=1e-12)  # correct to consider absolute tolerance only?

    def test_shear(self):
        # test whether shear is also applied (y is the axial direction)
        print(f"Working directory for test_shear: {self.parameters['settings']['working_directory']}")

        # create solver
        solver = create_instance(self.parameters)
        interface_input = solver.get_interface_input()

        # define a non-zero shear in y-direction
        local_shear = self.shear + np.array([0, 5, 0])

        # give value to variables
        pressure = interface_input.get_variable_data(self.mp_name_in, 'pressure')
        pressure[:] = self.p
        interface_input.set_variable_data(self.mp_name_in, 'pressure', pressure)
        traction = interface_input.get_variable_data(self.mp_name_in, 'traction')
        traction[:, :] = local_shear
        interface_input.set_variable_data(self.mp_name_in, 'traction', traction)

        # do 4 steps
        solver.initialize()
        for i in range(4):
            solver.initialize_solution_step()
            solver.solve_solution_step(interface_input)
            solver.finalize_solution_step()
        solver.finalize()

        # compare output, as shear input has changed these should be different
        output_shear = solver.get_interface_output()
        self.a5 = output_shear.get_variable_data(self.mp_name_out, 'displacement').copy()

        # normalize data and compare
        self.mean_disp_y_shear = np.mean(self.a5[:, 1])/self.a5.shape[0]

        print(f'Mean y-displacement without shear = {self.mean_disp_y_no_shear} m')
        print(f'Mean y-displacement with shear = {self.mean_disp_y_shear} m')

        self.assertNotAlmostEqual(self.mean_disp_y_no_shear - self.mean_disp_y_shear, 0., delta=1e-12)


class TestSolverWrapperAbaqus614Tube3D(TestSolverWrapperAbaqus614Tube2D):
    setup_case = True
    file_name = join(os.path.dirname(__file__), 'test_v614/tube3d/parameters.json')
    dir_tmp = join(os.path.realpath(os.path.dirname(__file__)), 'test_v614/tube3d')
    p = 1500
    shear = np.array([0, 0, 0])
    mp_name_in = 'WALLOUTSIDE_load_points'
    mp_name_out = 'WALLOUTSIDE_nodes'

    def test_shear(self):
        # test whether shear is also applied (x is the axial direction)
        print(f"Working directory for test_shear: {self.parameters['settings']['working_directory']}")

        # create solver
        solver = create_instance(self.parameters)
        interface_input = solver.get_interface_input()

        # define a non-zero shear in x-direction
        local_shear = self.shear + np.array([5, 0, 0])

        # give value to variables
        pressure = interface_input.get_variable_data(self.mp_name_in, 'pressure')
        pressure[:] = self.p
        interface_input.set_variable_data(self.mp_name_in, 'pressure', pressure)
        traction = interface_input.get_variable_data(self.mp_name_in, 'traction')
        traction[:, :] = local_shear
        interface_input.set_variable_data(self.mp_name_in, 'traction', traction)

        # do 4 steps
        solver.initialize()
        for i in range(4):
            solver.initialize_solution_step()
            solver.solve_solution_step(interface_input)
            solver.finalize_solution_step()
        solver.finalize()

        # compare output, as shear input has changed these should be different
        output_shear = solver.get_interface_output()
        self.a5 = output_shear.get_variable_data(self.mp_name_out, 'displacement').copy()

        # normalize data and compare
        self.mean_disp_x_shear = np.mean(self.a5[:, 0])/self.a5.shape[0]

        print(f'Mean x-displacement without shear = {self.mean_disp_x_no_shear} m')
        print(f'Mean x-displacement with shear = {self.mean_disp_x_shear} m')

        self.assertNotAlmostEqual(self.mean_disp_x_no_shear - self.mean_disp_x_shear, 0., delta=1e-12)


if __name__ == '__main__':
    unittest.main()
