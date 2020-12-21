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

            # Create the solver (__init__)
            AbaqusSolver0 = CreateInstance(par_solver_0)

            # give value to variables
            mp = AbaqusSolver0.model['BEAMINSIDEMOVING_load_points']
            for node in mp.Nodes:
                # domain extends from Y -0.025 to 0.025, default x-position is 0.005
                node.SetSolutionStepValue(pressure, 0, p)
                node.SetSolutionStepValue(traction, 0, [shear_x, shear_y, shear_z])

            AbaqusSolver0.Initialize()

            # Step 1, Coupling 1
            AbaqusSolver0.InitializeSolutionStep()
            output1_1 = AbaqusSolver0.SolveSolutionStep(AbaqusSolver0.GetInterfaceInput())
            # TODO: this could break down, depending on cwd check
            os.system('cp -r test_v614/tube2d/CSM/CSM_Time1.odb test_v614/tube2d/CSM/CSM_Time1_Iter1.odb')
            # Step 1, Coupling 2
            output1_2 = AbaqusSolver0.SolveSolutionStep(AbaqusSolver0.GetInterfaceInput()).deepcopy()
            AbaqusSolver0.FinalizeSolutionStep()

            # compare output, as input hasn't changed these should be the same
            a1 = output1_1.GetNumpyArray()
            a2 = output1_2.GetNumpyArray()

            # normalize data and compare
            mean = np.mean(a1)
            ref = np.abs(a1 - mean).max()

            a1n = (a1 - mean) / ref
            a2n = (a2 - mean) / ref

            for i in range(a1.size):
                self.assertAlmostEqual(a1n[i] - a2n[i], 0., delta=1e-12)

            # Step 2 and 3
            for i in range(2):
                AbaqusSolver0.InitializeSolutionStep()
                AbaqusSolver0.SolveSolutionStep(AbaqusSolver0.GetInterfaceInput())
                AbaqusSolver0.FinalizeSolutionStep()

            # Step 4
            AbaqusSolver0.InitializeSolutionStep()
            output_single_run = AbaqusSolver0.SolveSolutionStep(AbaqusSolver0.GetInterfaceInput()).deepcopy()
            AbaqusSolver0.FinalizeSolutionStep()
            AbaqusSolver0.Finalize()
            # TODO: this could break down, depending on cwd check
            os.system('cp test_v614/tube2d/CSM/CSM_Time4Surface0Output.dat'
                      ' test_v614/tube2d/CSM/CSM_Time4Surface0Output_Single.dat')

            # With restart
            # create solver which restarts at time step 2
            par_solver['settings'].SetInt('timestep_start', 2)
            AbaqusSolver1 = CreateInstance(par_solver)
            mp = AbaqusSolver1.model['BEAMINSIDEMOVING_load_points']
            for node in mp.Nodes:
                # Domain extends from Y -0.025 to 0.025, default x-position is 0.005
                node.SetSolutionStepValue(pressure, 0, p)
                node.SetSolutionStepValue(traction, 0, [shear_x, shear_y, shear_z])

            AbaqusSolver1.Initialize()

            # do step 3 and 4
            for i in range(2):
                AbaqusSolver1.InitializeSolutionStep()
                output_restart = AbaqusSolver1.SolveSolutionStep(AbaqusSolver1.GetInterfaceInput()).deepcopy()
                AbaqusSolver1.FinalizeSolutionStep()
            AbaqusSolver1.Finalize()

            # Compare output, as input hasn't changed these should be the same
            a1 = output_single_run.GetNumpyArray()
            a2 = output_restart.GetNumpyArray()

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
            AbaqusSolver2 = CreateInstance(par_solver)

            # give value to variables
            mp = AbaqusSolver2.model['BEAMINSIDEMOVING_load_points']
            for node in mp.Nodes:
                # Domain extends from Y -0.025 to 0.025, default x-position is 0.005
                node.SetSolutionStepValue(pressure, 0, p)
                node.SetSolutionStepValue(traction, 0, [shear_x, shear_y, shear_z])

            # do 4 steps
            AbaqusSolver2.Initialize()
            for i in range(4):
                AbaqusSolver2.InitializeSolutionStep()
                output_4cores = AbaqusSolver2.SolveSolutionStep(AbaqusSolver2.GetInterfaceInput()).deepcopy()
                AbaqusSolver2.FinalizeSolutionStep()
            AbaqusSolver2.Finalize()

            # Compare output, as input hasn't changed these should be the same
            # normalize data and compare
            a4 = output_4cores.GetNumpyArray()
            a4n = (a4 - mean)/ref

            for i in range(a1.size):
                self.assertAlmostEqual(a2n[i] - a4n[i], 0., delta=1e-12)
                self.assertAlmostEqual(a1n[i] - a4n[i], 0., delta=1e-12)

        # test whether shear is also applied
        if True:
            print_box('test whether shear is also applied')
            shear_y = 5

            # give value to variables
            mp = AbaqusSolver2.model['BEAMINSIDEMOVING_load_points']
            for node in mp.Nodes:
                # Domain extends from Y -0.025 to 0.025, default x-position is 0.005
                node.SetSolutionStepValue(pressure, 0, p)
                node.SetSolutionStepValue(traction, 0, [shear_x, shear_y, shear_z])

            # do 4 steps
            AbaqusSolver2.Initialize()
            for i in range(4):
                AbaqusSolver2.InitializeSolutionStep()
                output_shear = AbaqusSolver2.SolveSolutionStep(AbaqusSolver2.GetInterfaceInput()).deepcopy()
                AbaqusSolver2.FinalizeSolutionStep()
            AbaqusSolver2.Finalize()

            # Compare output, as shear input has changed these should be different
            # normalize data and compare
            a5 = output_shear.GetNumpyArray()

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
        # Axial direction is the x-direction

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
            p_setup_abaqus = subprocess.Popen(os.path.join(dir_tmp, 'setup_abaqus.sh'), cwd=dir_tmp, shell=True)
            p_setup_abaqus.wait()

        # test start and restart
        if True:
            # Create the solver (__init__)
            AbaqusSolver0 = CreateInstance(par_solver_0)

            # give value to variables
            mp = AbaqusSolver0.model['WALLOUTSIDE_load_points']
            for node in mp.Nodes:
                # Domain extends from Y -0.025 to 0.025, default x-position is 0.005
                node.SetSolutionStepValue(pressure, 0, p)
                node.SetSolutionStepValue(traction, 0, [shear_x, shear_y, shear_z])

            AbaqusSolver0.Initialize()

            # Step 1, Coupling 1
            AbaqusSolver0.InitializeSolutionStep()
            output1_1 = AbaqusSolver0.SolveSolutionStep(AbaqusSolver0.GetInterfaceInput())
            os.system('cp -r test_v614/tube3d/CSM/CSM_Time1.odb test_v614/tube3d/CSM/CSM_Time1_Iter1.odb')
            # Step 1, Coupling 2
            output1_2 = AbaqusSolver0.SolveSolutionStep(AbaqusSolver0.GetInterfaceInput()).deepcopy()
            AbaqusSolver0.FinalizeSolutionStep()

            # Compare output, as input hasn't changed these should be the same
            a1 = output1_1.GetNumpyArray()
            a2 = output1_2.GetNumpyArray()

            # normalize data and compare
            mean = np.mean(a1)
            ref = np.abs(a1 - mean).max()

            a1n = (a1 - mean) / ref
            a2n = (a2 - mean) / ref

            for i in range(a1.size):
                self.assertAlmostEqual(a1n[i] - a2n[i], 0., delta=1e-12)

            # Step 2 and 3
            for i in range(2):
                AbaqusSolver0.InitializeSolutionStep()
                AbaqusSolver0.SolveSolutionStep(AbaqusSolver0.GetInterfaceInput())
                AbaqusSolver0.FinalizeSolutionStep()
            # Step 4
            AbaqusSolver0.InitializeSolutionStep()
            output_single_run = AbaqusSolver0.SolveSolutionStep(AbaqusSolver0.GetInterfaceInput()).deepcopy()
            AbaqusSolver0.FinalizeSolutionStep()
            AbaqusSolver0.Finalize()

            os.system('cp test_v614/tube3d/CSM/CSM_Time4Surface0Output.dat '
                      'test_v614/tube3d/CSM/CSM_Time4Surface0Output_Single.dat')

            # With restart
            # create solver which restarts at time step 2
            par_solver['settings'].SetInt('timestep_start', 2)
            AbaqusSolver1 = CreateInstance(par_solver)

            # give value to variables
            mp = AbaqusSolver1.model['WALLOUTSIDE_load_points']
            for node in mp.Nodes:
                # Domain extends from Y -0.025 to 0.025, default x-position is 0.005
                node.SetSolutionStepValue(pressure, 0, p)
                node.SetSolutionStepValue(traction, 0, [shear_x, shear_y, shear_z])

            AbaqusSolver1.Initialize()

            # do step 3 and 4
            for i in range(2):
                AbaqusSolver1.InitializeSolutionStep()
                output_restart = AbaqusSolver1.SolveSolutionStep(AbaqusSolver1.GetInterfaceInput()).deepcopy()
                AbaqusSolver1.FinalizeSolutionStep()
            AbaqusSolver1.Finalize()

            # Compare output, as input hasn't changed these should be the same
            a1 = output_single_run.GetNumpyArray()
            a2 = output_restart.GetNumpyArray()

            # normalize data and compare
            mean = np.mean(a1)
            ref = np.abs(a1 - mean).max()

            a1n = (a1 - mean) / ref
            a2n = (a2 - mean) / ref

            for i in range(a1.size):
                self.assertAlmostEqual(a1n[i] - a2n[i], 0., delta=1e-12)

        # test whether using 4 CPUs gives the same results as using a single one
        if True:
            print_box('test whether using 4 CPUs gives the same results as using a single one.')

            # adapt Parameters, create solver
            par_solver = deepcopy(par_solver_0)
            par_solver['settings'].SetInt('cores', 4)
            AbaqusSolver2 = CreateInstance(par_solver)

            # give value to variables
            mp = AbaqusSolver2.model['WALLOUTSIDE_load_points']
            for node in mp.Nodes:
                # Domain extends from Y -0.025 to 0.025, default x-position is 0.005
                node.SetSolutionStepValue(pressure, 0, p)
                node.SetSolutionStepValue(traction, 0, [shear_x, shear_y, shear_z])

            # do 4 steps
            AbaqusSolver2.Initialize()
            for i in range(4):
                AbaqusSolver2.InitializeSolutionStep()
                output_4cores = AbaqusSolver2.SolveSolutionStep(AbaqusSolver2.GetInterfaceInput()).deepcopy()
                AbaqusSolver2.FinalizeSolutionStep()
            AbaqusSolver2.Finalize()

            # Compare output, as input hasn't changed these should be the same
            # normalize data and compare
            a4 = output_4cores.GetNumpyArray()
            a4n = (a4 - mean) / ref

            for i in range(a1.size):
                self.assertAlmostEqual(a2n[i] - a4n[i], 0., delta=1e-12)
                self.assertAlmostEqual(a1n[i] - a4n[i], 0., delta=1e-12)

        # test whether shear is also applied (x is the axial direction)
        if True:
            print_box('test whether shear is also applied')

            shear_x = 5
            mp = AbaqusSolver2.model['WALLOUTSIDE_load_points']

            # give value to variables
            for node in mp.Nodes:
                # Domain extends from Y -0.025 to 0.025, default x-position is 0.005
                node.SetSolutionStepValue(pressure, 0, p)
                node.SetSolutionStepValue(traction, 0, [shear_x, shear_y, shear_z])

            # do 4 steps
            AbaqusSolver2.Initialize()
            for i in range(4):
                AbaqusSolver2.InitializeSolutionStep()
                output_shear = AbaqusSolver2.SolveSolutionStep(AbaqusSolver2.GetInterfaceInput()).deepcopy()
                AbaqusSolver2.FinalizeSolutionStep()
            AbaqusSolver2.Finalize()

            # Compare output, as shear input has changed these should be different
            # normalize data and compare
            a5 = output_shear.GetNumpyArray()

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
