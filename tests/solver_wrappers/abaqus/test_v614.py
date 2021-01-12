from coconut.coupling_components.tools import create_instance

import unittest
import numpy as np
import os
from os.path import join
import subprocess
import json


def print_box(text):
    n = len(text)
    top = '\n┌─' + n * '─' + '─┐'
    mid = '\n│ ' + text + ' │'
    bottom = '\n└─' + n * '─' + '─┘'
    print(top + mid + bottom)


class TestSolverWrapperAbaqus614Tube2D(unittest.TestCase):
    setup_case = True

    """
    TODO: now, different test methods can not be run separately. I will take out the part from the test_restart method
    up to the actual restart and put that in the setUpClass method, as all other methods need that output to compare
    with later on. As such, the actual test methods are completely independent of each other."""

    @classmethod
    def setUpClass(cls):
        if cls.setup_case:
            dir_tmp = join(os.path.realpath(os.path.dirname(__file__)), 'test_v614/tube2d')
            p_setup_abaqus = subprocess.Popen(os.path.join(dir_tmp, 'setup_abaqus.sh'), cwd=dir_tmp, shell=True)
            p_setup_abaqus.wait()

    def setUp(self):
        file_name = join(os.path.dirname(__file__), 'test_v614/tube2d/parameters.json')
        with open(file_name, 'r') as parameter_file:
            self.parameters = json.load(parameter_file)
        self.p = 1500
        self.shear = np.array([0, 0, 0])
        self.mp_name_in = 'BEAMINSIDEMOVING_load_points'
        self.mp_name_out = 'BEAMINSIDEMOVING_nodes'

    def test_restart(self):
        # test if restart option works correctly

        # create the solver
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
        # step 1, coupling 1
        solver.initialize_solution_step()
        output1_1 = solver.solve_solution_step(interface_input)
        # step 1, coupling 2
        output1_2 = solver.solve_solution_step(interface_input)
        solver.finalize_solution_step()

        # compare output, as input hasn't changed these should be the same
        a1 = output1_1.get_variable_data(self.mp_name_out, 'displacement').copy()
        a2 = output1_2.get_variable_data(self.mp_name_out, 'displacement').copy()

        # normalize data and compare
        mean = np.mean(a1)
        ref = np.abs(a1 - mean).max()
        a1n = (a1 - mean) / ref
        a2n = (a2 - mean) / ref
        np.testing.assert_allclose(a1n, a2n, rtol=1e-12)

        # step 2 to 4
        for i in range(3):
            solver.initialize_solution_step()
            solver.solve_solution_step(interface_input)
            solver.finalize_solution_step()
        solver.finalize()

        # get data for solver without restart
        output_single_run = solver.get_interface_output()
        a1 = output_single_run.get_variable_data(self.mp_name_out, 'displacement').copy()

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
        a2 = output_restart.get_variable_data(self.mp_name_out, 'displacement').copy()
        self.mean_disp_y_no_shear = np.mean(a2[:, 1])/a2.shape[1]

        # normalize data and compare
        self.mean = np.mean(a1)
        self.ref = np.abs(a1 - mean).max()
        self.a1n = (a1 - self.mean) / self.ref
        self.a2n = (a2 - self.mean) / self.ref

        np.testing.assert_allclose(self.a1n, self.a2n, rtol=1e-12)

    def test_partitioning(self):
        # test whether using 4 CPUs gives the same results as using a single one

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
        # normalize data and compare
        output_4cores = solver.get_interface_output()
        a4 = output_4cores.get_variable_data(self.mp_name_out, 'displacement').copy()
        self.a4n = (a4 - self.mean) / self.ref

        np.testing.assert_allclose(self.a4n, self.a2n, rtol=1e-12)
        np.testing.assert_allclose(self.a4n, self.a1n, rtol=1e-12)

    def test_shear(self):
        # test whether shear is also applied (y is the axial direction)

        # create solver
        solver = create_instance(self.parameters)
        interface_input = solver.get_interface_input()

        # define a non-zero shear in y-direction
        self.shear[1] = 5

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

        # compare output, as shear input has changed these should be different
        output_shear = solver.get_interface_output
        a5 = output_shear.get_variable_data(self.mp_name_out, 'displacement').copy()

        # normalize data and compare
        self.mean_disp_y_shear = np.mean(a5[:, 1])/a5.shape[1]

        print(f'Mean y-displacement without shear = {self.mean_disp_y_no_shear} m')
        print(f'Mean y-displacement with shear = {self.mean_disp_y_shear} m')

        self.assertNotAlmostEqual(self.mean_disp_y_no_shear - self.mean_disp_y_shear, 0., delta=1e-12)


# TODO: put in separate class and adapt
# class TestSolverWrapperAbaqus614Tube3D(unittest.TestCase):
#     print_box('started tests for Abaqus Tube3D')
#     # axial direction is the x-direction
#
#     parameter_file_name = os.path.join(os.path.dirname(__file__), 'test_v614/tube3d', 'parameters.json')
#
#     with open(parameter_file_name, 'r') as parameter_file:
#         parameters = data_structure.Parameters(parameter_file.read())
#     par_solver_0 = parameters['solver_wrappers'][0]
#
#     # if running from this folder
#     if os.getcwd() == os.path.realpath(os.path.dirname(__file__)):
#         par_solver_0['settings'].SetString('working_directory', 'test_v614/tube3d/CSM')
#
#     par_solver = deepcopy(par_solver_0)
#
#     # "global" definitions
#     pressure = vars(data_structure)['PRESSURE']
#     traction = vars(data_structure)['TRACTION']
#
#     p = 1500
#     shear_x = 0
#     shear_y = 0
#     shear_z = 0
#
#     # setup Abaqus case
#     if True:
#         print_box('setup Abaqus case')
#         dir_tmp = os.path.join(os.path.realpath(os.path.dirname(__file__)), f'test_v614/tube3d')
#         print(f'dir_tmp = {dir_tmp}')
#         p_setup_abaqus = subprocess.Popen(os.path.join(dir_tmp, 'setup_abaqus.sh'), cwd=dir_tmp, shell=True)
#         p_setup_abaqus.wait()
#
#     # test start and restart
#     if True:
#         print_box('test start and restart')
#
#         # create the solver (__init__)
#         solver = CreateInstance(par_solver)
#
#         # give value to variables
#         mp = solver.model['WALLOUTSIDE_load_points']
#         for node in mp.Nodes:
#             # domain extends from Y -0.025 to 0.025, default x-position is 0.005
#             node.SetSolutionStepValue(pressure, 0, p)
#             node.SetSolutionStepValue(traction, 0, [shear_x, shear_y, shear_z])
#
#         solver.Initialize()
#
#         # step 1, coupling 1
#         solver.InitializeSolutionStep()
#         output1_1 = solver.SolveSolutionStep(solver.GetInterfaceInput()).deepcopy()
#         # step 1, coupling 2
#         output1_2 = solver.SolveSolutionStep(solver.GetInterfaceInput()).deepcopy()
#         solver.FinalizeSolutionStep()
#
#         # compare output, as input hasn't changed these should be the same
#         a1 = output1_1.GetNumpyArray().copy()
#         a2 = output1_2.GetNumpyArray().copy()
#
#         # normalize data and compare
#         mean = np.mean(a1)
#         ref = np.abs(a1 - mean).max()
#
#         a1n = (a1 - mean) / ref
#         a2n = (a2 - mean) / ref
#
#         for i in range(a1.size):
#             self.assertAlmostEqual(a1n[i] - a2n[i], 0., delta=1e-12)
#
#         # step 2 to 4
#         for i in range(3):
#             solver.InitializeSolutionStep()
#             solver.SolveSolutionStep(solver.GetInterfaceInput())
#             solver.FinalizeSolutionStep()
#         solver.Finalize()
#
#         # get data for solver without restart
#         output_single_run = solver.GetInterfaceOutput().deepcopy()
#         a1 = output_single_run.GetNumpyArray().copy()
#
#         # create solver which restarts at time step 2
#         par_solver['settings'].SetInt('timestep_start', 2)
#         solver = CreateInstance(par_solver)
#
#         mp = solver.model['WALLOUTSIDE_load_points']
#
#         # give value to variables
#         for node in mp.Nodes:
#             # domain extends from Y -0.025 to 0.025, default x-position is 0.005
#             node.SetSolutionStepValue(pressure, 0, p)
#             node.SetSolutionStepValue(traction, 0, [shear_x, shear_y, shear_z])
#
#         solver.Initialize()
#
#         # do step 3 and 4
#         for i in range(2):
#             solver.InitializeSolutionStep()
#             solver.SolveSolutionStep(solver.GetInterfaceInput())
#             solver.FinalizeSolutionStep()
#         solver.Finalize()
#
#         # compare output, as input hasn't changed these should be the same
#         # get data for solver with restart
#         output_restart = solver.GetInterfaceOutput().deepcopy()
#         a2 = output_restart.GetNumpyArray().copy()
#
#         # normalize data and compare
#         mean = np.mean(a1)
#         ref = np.abs(a1 - mean).max()
#
#         a1n = (a1 - mean) / ref
#         a2n = (a2 - mean) / ref
#
#         for i in range(a1.size):
#             self.assertAlmostEqual(a1n[i] - a2n[i], 0., delta=1e-12)
#
#     # test whether using 4 CPUs gives the same results as using a single one
#     if True:
#         print_box('test whether using 4 CPUs gives the same results as using a single one')
#
#         # adapt Parameters, create solver
#         par_solver = deepcopy(par_solver_0)
#         par_solver['settings'].SetInt('cores', 4)
#         solver = CreateInstance(par_solver)
#
#         # give value to variables
#         mp = solver.model['WALLOUTSIDE_load_points']
#         for node in mp.Nodes:
#             # domain extends from Y -0.025 to 0.025, default x-position is 0.005
#             node.SetSolutionStepValue(pressure, 0, p)
#             node.SetSolutionStepValue(traction, 0, [shear_x, shear_y, shear_z])
#
#         # do 4 steps
#         solver.Initialize()
#         for i in range(4):
#             solver.InitializeSolutionStep()
#             solver.SolveSolutionStep(solver.GetInterfaceInput())
#             solver.FinalizeSolutionStep()
#         solver.Finalize()
#
#         # compare output, as input hasn't changed these should be the same
#         # normalize data and compare
#         output_4cores = solver.GetInterfaceOutput().deepcopy()
#         a4 = output_4cores.GetNumpyArray().copy()
#         a4n = (a4 - mean) / ref
#
#         for i in range(a1.size):
#             self.assertAlmostEqual(a2n[i] - a4n[i], 0., delta=1e-12)
#             self.assertAlmostEqual(a1n[i] - a4n[i], 0., delta=1e-12)
#
#     # test whether shear is also applied (x is the axial direction)
#     if True:
#         print_box('test whether shear is also applied')
#
#         # create solver
#         par_solver = deepcopy(par_solver_0)
#         solver = CreateInstance(par_solver)
#
#         shear_x = 5
#
#         mp = solver.model['WALLOUTSIDE_load_points']
#
#         # give value to variables
#         for node in mp.Nodes:
#             # domain extends from Y -0.025 to 0.025, default x-position is 0.005
#             node.SetSolutionStepValue(pressure, 0, p)
#             node.SetSolutionStepValue(traction, 0, [shear_x, shear_y, shear_z])
#
#         # do 4 steps
#         solver.Initialize()
#         for i in range(4):
#             solver.InitializeSolutionStep()
#             solver.SolveSolutionStep(solver.GetInterfaceInput())
#             solver.FinalizeSolutionStep()
#         solver.Finalize()
#
#         # compare output, as shear input has changed these should be different
#         output_shear = solver.GetInterfaceOutput().deepcopy()
#         a5 = output_shear.GetNumpyArray().copy()
#
#         # normalize data and compare
#         mean_disp_x_no_shear = 0
#         mean_disp_x_shear = 0
#
#         for i in range(0, a1.size, 3):
#             mean_disp_x_no_shear += a2[i]
#             mean_disp_x_shear += a5[i]
#
#         mean_disp_x_no_shear /= (a1.size / 3)
#         mean_disp_x_shear /= (a1.size / 3)
#
#         print(f'Mean x-displacement without shear = {mean_disp_x_no_shear} m')
#         print(f'Mean x-displacement with shear = {mean_disp_x_shear} m')
#
#         self.assertNotAlmostEqual(mean_disp_x_no_shear - mean_disp_x_shear, 0., delta=1e-12)


if __name__ == '__main__':
    unittest.main()
