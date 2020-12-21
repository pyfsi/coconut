from coconut import data_structure
from coconut.data_structure import KratosUnittest
from coconut.coupling_components.tools import CreateInstance

import numpy as np
from copy import deepcopy
import os
import subprocess


def print_box(text):
    n = len(text)
    top = '\n┌─' + n * '─' + '─┐'
    mid = '\n│ ' + text + ' │'
    bottom = '\n└─' + n * '─' + '─┘'
    print(top + mid + bottom)


class TestSolverWrapperAbaqus614(KratosUnittest.TestCase):
    def test_solver_wrapper_abaqus_614(self):
        self.test_solver_wrapper_abaqus_614_tube2d()
        self.test_solver_wrapper_abaqus_614_tube3d()

    def test_solver_wrapper_abaqus_614_tube2d(self):
        print_box('started tests for Abaqus Tube2D')
        # axial direction is the y-direction

        parameter_file_name = os.path.join(os.path.dirname(__file__), 'test_v614/tube2d', 'test_solver_wrapper.json')

        with open(parameter_file_name, 'r') as parameter_file:
            parameters = data_structure.Parameters(parameter_file.read())
        par_solver_0 = parameters['solver_wrappers'][0]

        # if running from this folder
        if os.getcwd() == os.path.realpath(os.path.dirname(__file__)):
            par_solver_0['settings'].SetString('working_directory', 'test_v614/tube2d/CSM')

        par_solver = deepcopy(par_solver_0)

        # "global" definitions
        pressure = vars(data_structure)['PRESSURE']
        traction = vars(data_structure)['TRACTION']

        p = 1500
        shear_x = 0
        shear_y = 0
        shear_z = 0

        # setup Abaqus case
        if True:
            print_box('setup Abaqus case')
            dir_tmp = os.path.join(os.path.realpath(os.path.dirname(__file__)), f'test_v614/tube2d')
            print(f'dir_tmp = {dir_tmp}')
            p_setup_abaqus = subprocess.Popen(os.path.join(dir_tmp, 'setup_abaqus.sh'), cwd=dir_tmp, shell=True)
            p_setup_abaqus.wait()

        # test start and restart
        if True:
            print_box('test start and restart')

            # create the solver (__init__)
            solver = CreateInstance(par_solver)

            # give value to variables
            mp = solver.model['BEAMINSIDEMOVING_load_points']
            for node in mp.Nodes:
                # domain extends from Y -0.025 to 0.025, default x-position is 0.005
                node.SetSolutionStepValue(pressure, 0, p)
                node.SetSolutionStepValue(traction, 0, [shear_x, shear_y, shear_z])

            solver.Initialize()

            # step 1, coupling 1
            solver.InitializeSolutionStep()
            output1_1 = solver.SolveSolutionStep(solver.GetInterfaceInput()).deepcopy()
            # step 1, coupling 2
            output1_2 = solver.SolveSolutionStep(solver.GetInterfaceInput()).deepcopy()
            solver.FinalizeSolutionStep()

            # compare output, as input hasn't changed these should be the same
            a1 = output1_1.GetNumpyArray().copy()
            a2 = output1_2.GetNumpyArray().copy()

            # normalize data and compare
            mean = np.mean(a1)
            ref = np.abs(a1 - mean).max()

            a1n = (a1 - mean) / ref
            a2n = (a2 - mean) / ref

            for i in range(a1.size):
                self.assertAlmostEqual(a1n[i] - a2n[i], 0., delta=1e-12)

            # step 2 to 4
            for i in range(3):
                solver.InitializeSolutionStep()
                solver.SolveSolutionStep(solver.GetInterfaceInput())
                solver.FinalizeSolutionStep()
            solver.Finalize()

            # get data for solver without restart
            output_single_run = solver.SolveSolutionStep(solver.GetInterfaceInput()).deepcopy()
            a1 = output_single_run.GetNumpyArray().copy()

            # create solver which restarts at time step 2
            par_solver['settings'].SetInt('timestep_start', 2)
            solver = CreateInstance(par_solver)

            mp = solver.model['BEAMINSIDEMOVING_load_points']

            # give value to variables
            for node in mp.Nodes:
                # domain extends from Y -0.025 to 0.025, default x-position is 0.005
                node.SetSolutionStepValue(pressure, 0, p)
                node.SetSolutionStepValue(traction, 0, [shear_x, shear_y, shear_z])

            solver.Initialize()

            # do step 3 and 4
            for i in range(2):
                solver.InitializeSolutionStep()
                solver.SolveSolutionStep(solver.GetInterfaceInput())
                solver.FinalizeSolutionStep()
            solver.Finalize()

            # compare output, as input hasn't changed these should be the same
            # get data for solver with restart
            output_restart = solver.SolveSolutionStep(solver.GetInterfaceInput()).deepcopy()
            a2 = output_restart.GetNumpyArray().copy()

            # normalize data and compare
            mean = np.mean(a1)
            ref = np.abs(a1 - mean).max()

            a1n = (a1 - mean) / ref
            a2n = (a2 - mean) / ref

            for i in range(a1.size):
                self.assertAlmostEqual(a1n[i] - a2n[i], 0., delta=1e-12)

        # test whether using 4 CPUs gives the same results as using a single one
        if True:
            print_box('test whether using 4 CPUs gives the same results as using a single one')

            # adapt Parameters, create solver
            par_solver = deepcopy(par_solver_0)
            par_solver['settings'].SetInt('cores', 4)
            solver = CreateInstance(par_solver)

            # give value to variables
            mp = solver.model['BEAMINSIDEMOVING_load_points']
            for node in mp.Nodes:
                # domain extends from Y -0.025 to 0.025, default x-position is 0.005
                node.SetSolutionStepValue(pressure, 0, p)
                node.SetSolutionStepValue(traction, 0, [shear_x, shear_y, shear_z])

            # do 4 steps
            solver.Initialize()
            for i in range(4):
                solver.InitializeSolutionStep()
                solver.SolveSolutionStep(solver.GetInterfaceInput())
                solver.FinalizeSolutionStep()
            solver.Finalize()

            # compare output, as input hasn't changed these should be the same
            # normalize data and compare
            output_4cores = solver.SolveSolutionStep(solver.GetInterfaceInput()).deepcopy()
            a4 = output_4cores.GetNumpyArray().copy()
            a4n = (a4 - mean) / ref

            for i in range(a1.size):
                self.assertAlmostEqual(a2n[i] - a4n[i], 0., delta=1e-12)
                self.assertAlmostEqual(a1n[i] - a4n[i], 0., delta=1e-12)

        # test whether shear is also applied (y is the axial direction)
        if True:
            print_box('test whether shear is also applied')

            # create solver
            par_solver = deepcopy(par_solver_0)
            solver = CreateInstance(par_solver)

            shear_y = 5

            mp = solver.model['BEAMINSIDEMOVING_load_points']

            # give value to variables
            for node in mp.Nodes:
                # domain extends from Y -0.025 to 0.025, default x-position is 0.005
                node.SetSolutionStepValue(pressure, 0, p)
                node.SetSolutionStepValue(traction, 0, [shear_x, shear_y, shear_z])

            # do 4 steps
            solver.Initialize()
            for i in range(4):
                solver.InitializeSolutionStep()
                solver.SolveSolutionStep(solver.GetInterfaceInput())
                solver.FinalizeSolutionStep()
            solver.Finalize()

            # compare output, as shear input has changed these should be different
            output_shear = solver.SolveSolutionStep(solver.GetInterfaceInput()).deepcopy()
            a5 = output_shear.GetNumpyArray().copy()

            # normalize data and compare
            mean_disp_y_no_shear = 0
            mean_disp_y_shear = 0

            for i in range(1, a1.size, 3):
                mean_disp_y_no_shear += a2[i]
                mean_disp_y_shear += a5[i]

            mean_disp_y_no_shear /= (a1.size/3)
            mean_disp_y_shear /= (a1.size/3)

            print(f'Mean y-displacement without shear = {mean_disp_y_no_shear} m')
            print(f'Mean y-displacement with shear = {mean_disp_y_shear} m')

            self.assertNotAlmostEqual(mean_disp_y_no_shear - mean_disp_y_shear, 0., delta=1e-12)

    def test_solver_wrapper_abaqus_614_tube3d(self):
        print_box('started tests for Abaqus Tube3D')
        # axial direction is the x-direction

        parameter_file_name = os.path.join(os.path.dirname(__file__), 'test_v614/tube3d', 'test_solver_wrapper.json')

        with open(parameter_file_name, 'r') as parameter_file:
            parameters = data_structure.Parameters(parameter_file.read())
        par_solver_0 = parameters['solver_wrappers'][0]

        # if running from this folder
        if os.getcwd() == os.path.realpath(os.path.dirname(__file__)):
            par_solver_0['settings'].SetString('working_directory', 'test_v614/tube3d/CSM')

        par_solver = deepcopy(par_solver_0)

        # "global" definitions
        pressure = vars(data_structure)['PRESSURE']
        traction = vars(data_structure)['TRACTION']

        p = 1500
        shear_x = 0
        shear_y = 0
        shear_z = 0

        # setup Abaqus case
        if True:
            print_box('setup Abaqus case')
            dir_tmp = os.path.join(os.path.realpath(os.path.dirname(__file__)), f'test_v614/tube3d')
            print(f'dir_tmp = {dir_tmp}')
            p_setup_abaqus = subprocess.Popen(os.path.join(dir_tmp, 'setup_abaqus.sh'), cwd=dir_tmp, shell=True)
            p_setup_abaqus.wait()

        # test start and restart
        if True:
            print_box('test start and restart')

            # create the solver (__init__)
            solver = CreateInstance(par_solver)

            # give value to variables
            mp = solver.model['WALLOUTSIDE_load_points']
            for node in mp.Nodes:
                # domain extends from Y -0.025 to 0.025, default x-position is 0.005
                node.SetSolutionStepValue(pressure, 0, p)
                node.SetSolutionStepValue(traction, 0, [shear_x, shear_y, shear_z])

            solver.Initialize()

            # step 1, coupling 1
            solver.InitializeSolutionStep()
            output1_1 = solver.SolveSolutionStep(solver.GetInterfaceInput()).deepcopy()
            # step 1, coupling 2
            output1_2 = solver.SolveSolutionStep(solver.GetInterfaceInput()).deepcopy()
            solver.FinalizeSolutionStep()

            # compare output, as input hasn't changed these should be the same
            a1 = output1_1.GetNumpyArray().copy()
            a2 = output1_2.GetNumpyArray().copy()

            # normalize data and compare
            mean = np.mean(a1)
            ref = np.abs(a1 - mean).max()

            a1n = (a1 - mean) / ref
            a2n = (a2 - mean) / ref

            for i in range(a1.size):
                self.assertAlmostEqual(a1n[i] - a2n[i], 0., delta=1e-12)

            # step 2 to 4
            for i in range(3):
                solver.InitializeSolutionStep()
                solver.SolveSolutionStep(solver.GetInterfaceInput())
                solver.FinalizeSolutionStep()
            solver.Finalize()

            # get data for solver without restart
            output_single_run = solver.SolveSolutionStep(solver.GetInterfaceInput()).deepcopy()
            a1 = output_single_run.GetNumpyArray().copy()

            # create solver which restarts at time step 2
            par_solver['settings'].SetInt('timestep_start', 2)
            solver = CreateInstance(par_solver)

            mp = solver.model['WALLOUTSIDE_load_points']

            # give value to variables
            for node in mp.Nodes:
                # domain extends from Y -0.025 to 0.025, default x-position is 0.005
                node.SetSolutionStepValue(pressure, 0, p)
                node.SetSolutionStepValue(traction, 0, [shear_x, shear_y, shear_z])

            solver.Initialize()

            # do step 3 and 4
            for i in range(2):
                solver.InitializeSolutionStep()
                solver.SolveSolutionStep(solver.GetInterfaceInput())
                solver.FinalizeSolutionStep()
            solver.Finalize()

            # compare output, as input hasn't changed these should be the same
            # get data for solver with restart
            output_restart = solver.SolveSolutionStep(solver.GetInterfaceInput()).deepcopy()
            a2 = output_restart.GetNumpyArray().copy()

            # normalize data and compare
            mean = np.mean(a1)
            ref = np.abs(a1 - mean).max()

            a1n = (a1 - mean) / ref
            a2n = (a2 - mean) / ref

            for i in range(a1.size):
                self.assertAlmostEqual(a1n[i] - a2n[i], 0., delta=1e-12)

        # test whether using 4 CPUs gives the same results as using a single one
        if True:
            print_box('test whether using 4 CPUs gives the same results as using a single one')

            # adapt Parameters, create solver
            par_solver = deepcopy(par_solver_0)
            par_solver['settings'].SetInt('cores', 4)
            solver = CreateInstance(par_solver)

            # give value to variables
            mp = solver.model['WALLOUTSIDE_load_points']
            for node in mp.Nodes:
                # domain extends from Y -0.025 to 0.025, default x-position is 0.005
                node.SetSolutionStepValue(pressure, 0, p)
                node.SetSolutionStepValue(traction, 0, [shear_x, shear_y, shear_z])

            # do 4 steps
            solver.Initialize()
            for i in range(4):
                solver.InitializeSolutionStep()
                solver.SolveSolutionStep(solver.GetInterfaceInput())
                solver.FinalizeSolutionStep()
            solver.Finalize()

            # compare output, as input hasn't changed these should be the same
            # normalize data and compare
            output_4cores = solver.SolveSolutionStep(solver.GetInterfaceInput()).deepcopy()
            a4 = output_4cores.GetNumpyArray().copy()
            a4n = (a4 - mean) / ref

            for i in range(a1.size):
                self.assertAlmostEqual(a2n[i] - a4n[i], 0., delta=1e-12)
                self.assertAlmostEqual(a1n[i] - a4n[i], 0., delta=1e-12)

        # test whether shear is also applied (x is the axial direction)
        if True:
            print_box('test whether shear is also applied')

            # create solver
            par_solver = deepcopy(par_solver_0)
            solver = CreateInstance(par_solver)

            shear_x = 5

            mp = solver.model['WALLOUTSIDE_load_points']

            # give value to variables
            for node in mp.Nodes:
                # domain extends from Y -0.025 to 0.025, default x-position is 0.005
                node.SetSolutionStepValue(pressure, 0, p)
                node.SetSolutionStepValue(traction, 0, [shear_x, shear_y, shear_z])

            # do 4 steps
            solver.Initialize()
            for i in range(4):
                solver.InitializeSolutionStep()
                solver.SolveSolutionStep(solver.GetInterfaceInput())
                solver.FinalizeSolutionStep()
            solver.Finalize()

            # compare output, as shear input has changed these should be different
            output_shear = solver.SolveSolutionStep(solver.GetInterfaceInput()).deepcopy()
            a5 = output_shear.GetNumpyArray().copy()

            # normalize data and compare
            mean_disp_x_no_shear = 0
            mean_disp_x_shear = 0

            for i in range(0, a1.size, 3):
                mean_disp_x_no_shear += a2[i]
                mean_disp_x_shear += a5[i]

            mean_disp_x_no_shear /= (a1.size / 3)
            mean_disp_x_shear /= (a1.size / 3)

            print(f'Mean x-displacement without shear = {mean_disp_x_no_shear} m')
            print(f'Mean x-displacement with shear = {mean_disp_x_shear} m')

            self.assertNotAlmostEqual(mean_disp_x_no_shear - mean_disp_x_shear, 0., delta=1e-12)


if __name__ == '__main__':
    KratosUnittest.main()
