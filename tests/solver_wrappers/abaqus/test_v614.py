from coconut.coupling_components.tools import create_instance

import unittest
import numpy as np
import os
from os.path import join
import subprocess
import json


class TestSolverWrapperAbaqus614Tube2D(unittest.TestCase):
    setup_case = True
    dimension = 2
    p = 1500
    shear = np.array([0, 0, 0])
    shear_dir = 1   # y-direction is axial direction
    mp_name_in = 'BEAMINSIDEMOVING_load_points'
    mp_name_out = 'BEAMINSIDEMOVING_nodes'

    @classmethod
    def setUpClass(cls):
        if cls.setup_case:
            dir_tmp = join(os.path.realpath(os.path.dirname(__file__)), f'test_v614/tube{int(cls.dimension)}d')
            p_setup_abaqus = subprocess.Popen(os.path.join(dir_tmp, 'setup_abaqus.sh'), cwd=dir_tmp, shell=True)
            p_setup_abaqus.wait()

        # perform reference calculation
        file_name = join(os.path.dirname(__file__), f'test_v614/tube{int(cls.dimension)}d/parameters.json')
        with open(file_name, 'r') as parameter_file:
            parameters = json.load(parameter_file)

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
        a1 = output1_1.get_variable_data(cls.mp_name_out, 'displacement')
        a2 = output1_2.get_variable_data(cls.mp_name_out, 'displacement')

        # compare
        np.testing.assert_allclose(a2, a1, rtol=1e-15)

        # step 2 to 4
        for i in range(3):
            solver.initialize_solution_step()
            solver.solve_solution_step(interface_input)
            solver.finalize_solution_step()
        solver.finalize()

        # get data for solver without restart
        output_single_run = solver.get_interface_output()
        cls.a1 = output_single_run.get_variable_data(cls.mp_name_out, 'displacement')
        print(f"Max disp a1: {np.max(np.abs(cls.a1), axis=0)}")

    def setUp(self):
        file_name = join(os.path.dirname(__file__), f'test_v614/tube{int(self.dimension)}d/parameters.json')
        with open(file_name, 'r') as parameter_file:
            self.parameters = json.load(parameter_file)

    def test_restart(self):
        """
        Test whether restarting at time step 2 and simulating 2 time steps yields the same displacement as when the
        simulation is ran from time step 0 until time step 4. A constant pressure is applied and no shear, on the tube
        examples. For the test the relative difference between displacements is checked, but it is required to also use
        a small absolute tolerance, otherwise the test will fail in the symmetry planes (i.e. whenever one of the
        original coordinates is 0), because those have a near-zero displacement after applying a uniform pressure an no
        shear. This near-zero value is a noise and has large relative differences between runs, but a very low absolute
        value, they are successfully filtered with an absolute tolerance.
        """

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

        # do step 3 and 4
        solver.initialize()
        for i in range(2):
            solver.initialize_solution_step()
            solver.solve_solution_step(interface_input)
            solver.finalize_solution_step()
        solver.finalize()

        # compare output, as input hasn't changed these should be the same
        # get data for solver with restart
        output_restart = solver.get_interface_output()
        self.a3 = output_restart.get_variable_data(self.mp_name_out, 'displacement')
        print(f"\nMax disp a3: {np.max(np.abs(self.a3), axis=0)}")
        print(f"Max diff between a1 and a3: {np.abs(self.a1 - self.a3).max(axis=0)}")

        if self.dimension == 2:
            np.testing.assert_array_equal(self.a3[:, 2], self.a3[:, 2] * 0.0)
            np.testing.assert_allclose(self.a3[:, :2], self.a1[:, :2], rtol=1e-10, atol=1e-17)
        else:
            np.testing.assert_allclose(self.a3, self.a1, rtol=1e-10, atol=1e-17)

    def test_partitioning(self):
        """
        Test whether using 4 CPUs yields the same results as using a single one. A constant pressure is applied and no
        shear, on the tube examples. For the test the relative difference between displacements is checked, but it is
        required to also use a small absolute tolerance, otherwise the test will fail in the symmetry planes (i.e.
        whenever one of the original coordinates is 0), because those have a near-zero displacement after applying a
        uniform pressure an no shear. This near-zero value is a noise and has large relative differences between runs,
        but a very low absolute value, they are successfully filtered with an absolute tolerance.
        """

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
        self.a4 = output_4cores.get_variable_data(self.mp_name_out, 'displacement')
        print(f"\nMax disp a4: {np.max(np.abs(self.a4), axis=0)}")
        print(f"Max diff between a1 and a4: {np.abs(self.a1 - self.a4).max(axis=0)}")

        if self.dimension == 2:
            np.testing.assert_array_equal(self.a4[:, 2], self.a4[:, 2] * 0.0)
            np.testing.assert_allclose(self.a4[:, :2], self.a1[:, :2], rtol=1e-10, atol=1e-17)
        else:
            np.testing.assert_allclose(self.a4, self.a1, rtol=1e-10, atol=1e-17)

    def test_shear(self):
        # test whether shear is also applied

        # create solver
        solver = create_instance(self.parameters)
        interface_input = solver.get_interface_input()

        # apply non-zero shear in axial_dir (y for 2D, x for 3D)
        local_shear = self.shear.copy()
        local_shear[self.axial_dir] = 5

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
        self.a5 = output_shear.get_variable_data(self.mp_name_out, 'displacement')
        self.mean_disp_shear = np.mean(np.abs(self.a5[:, self.shear_dir]))
        self.mean_disp_no_shear = np.mean(np.abs(self.a1[:, self.shear_dir]))
        print(f'Mean displacement in axial direction without shear = {self.mean_disp_no_shear} m')
        print(f'Mean displacement in axial direction with shear = {self.mean_disp_shear} m')
        self.assertNotAlmostEqual(self.mean_disp_no_shear - self.mean_disp_shear, 0., delta=1e-12)


class TestSolverWrapperAbaqus614Tube3D(TestSolverWrapperAbaqus614Tube2D):
    setup_case = True
    dimension = 3
    shear_dir = 0   # x-direction is axial direction
    mp_name_in = 'WALLOUTSIDE_load_points'
    mp_name_out = 'WALLOUTSIDE_nodes'


if __name__ == '__main__':
    unittest.main()
