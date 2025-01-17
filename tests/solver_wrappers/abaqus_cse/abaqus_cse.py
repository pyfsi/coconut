from coconut.tools import create_instance, rm_timed

import unittest
import numpy as np
import os
from os.path import join
import json
import shutil


class TestSolverWrapperAbaqusCSETube2D(unittest.TestCase):
    version = None  # Abaqus version, e.g. 2022, set in subclass
    setup_case = True
    dimension = 2
    axial_dir = 1  # y-direction is axial direction
    radial_dirs = [0]
    mp_name_in = 'BEAMINSIDEMOVING_load_points'
    mp_name_out = 'BEAMINSIDEMOVING_nodes'

    @classmethod
    def setUpClass(cls):
        dir_name = os.path.realpath(os.path.dirname(__file__))
        cls.file_name = join(dir_name, f'test_v{cls.version}/tube{cls.dimension}d/parameters.json')
        cls.working_dir = join(dir_name, f'test_v{cls.version}/tube{cls.dimension}d/CSM')
        cls.setup_dir = join(dir_name, f'test_v{cls.version}/tube{cls.dimension}d/setup_abaqus')

        # perform reference calculation
        with open(cls.file_name) as parameter_file:
            parameters = json.load(parameter_file)

        parameters['settings']['working_directory'] = os.path.relpath(cls.working_dir)  # set working directory

        if cls.setup_case:
            shutil.rmtree(cls.working_dir, ignore_errors=True)
            shutil.copytree(cls.setup_dir, cls.working_dir)

        # create the solver
        solver = create_instance(parameters)
        solver.initialize()
        interface_input = solver.get_interface_input()
        model_part = interface_input.get_model_part(cls.mp_name_in)
        coords = [model_part.x0, model_part.y0, model_part.z0]

        # give value to variables
        pressure = cls.get_p(coords[cls.axial_dir]).reshape(-1, 1)
        interface_input.set_variable_data(cls.mp_name_in, 'pressure', pressure)
        traction = np.zeros((model_part.size, 3))
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
            rm_timed(cls.working_dir)

    def setUp(self):
        with open(self.file_name) as parameter_file:
            self.parameters = json.load(parameter_file)
        self.parameters['settings']['working_directory'] = os.path.relpath(self.working_dir)  # set working directory

    @staticmethod
    def get_p(x):
        return 1500 * np.cos(2 * np.pi / 0.05 * x)

    @staticmethod
    def get_shear(x, axial_dir):
        shear = np.zeros((x.shape[0], 3))
        shear[:, axial_dir] = (x + 0.025) / 0.05 * 10
        return shear

    def test_repeat_iteration(self):
        """
        Test whether repeating the same iteration yields the same results.
        """
        np.testing.assert_allclose(self.a2_1, self.a1_1, rtol=1e-15)

    @unittest.skip('Not yet implemented, needs adjustment')
    def test_restart(self):
        """
        Test whether restarting at time step 2 and simulating 2 time steps yields the same displacement as when the
        simulation is run from time step 0 until time step 4. A constant pressure is applied and no shear, on the tube
        examples. For the test the relative difference between displacements is checked, but it is required to also use
        a small absolute tolerance, otherwise the test will fail in the symmetry planes (i.e. whenever one of the
        original coordinates is 0), because those have a near-zero displacement after applying a uniform pressure and no
        shear. This near-zero value is a noise and has large relative differences between runs, but a very low absolute
        value, they are successfully filtered with an absolute tolerance.
        """

        # create solver which restarts at time step 2
        self.parameters['settings']['timestep_start'] = 2
        solver = create_instance(self.parameters)
        solver.initialize()
        interface_input = solver.get_interface_input()
        model_part = interface_input.get_model_part(self.mp_name_in)
        coords = [model_part.x0, model_part.y0, model_part.z0]

        # give value to variables
        pressure = self.get_p(coords[self.axial_dir]).reshape(-1, 1)
        interface_input.set_variable_data(self.mp_name_in, 'pressure', pressure)
        traction = np.zeros((model_part.size, 3))
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
        indices = sorted([self.axial_dir] + self.radial_dirs)  # columns that contain non-zero data
        a3_extra = np.delete(self.a3, indices, axis=1)  # remove columns containing data
        np.testing.assert_array_equal(a3_extra, a3_extra * 0.0)  # if a column remains it should be all zeroes
        np.testing.assert_allclose(self.a3[:, indices], self.a1[:, indices], rtol=1e-10, atol=1e-17)  # non-zero columns

    def test_partitioning(self):
        """
        Test whether using 4 CPUs yields the same results as using a single one. A constant pressure is applied and no
        shear, on the tube examples. For the test the relative difference between displacements is checked, but it is
        required to also use a small absolute tolerance, otherwise the test will fail in the symmetry planes (i.e.
        whenever one of the original coordinates is 0), because those have a near-zero displacement after applying a
        uniform pressure and no shear. This near-zero value is a noise and has large relative differences between runs,
        but a very low absolute value, they are successfully filtered with an absolute tolerance.
        """

        # adapt Parameters, create solver
        self.parameters['settings']['cores'] = 4
        solver = create_instance(self.parameters)
        solver.initialize()
        interface_input = solver.get_interface_input()
        model_part = interface_input.get_model_part(self.mp_name_in)
        coords = [model_part.x0, model_part.y0, model_part.z0]

        # give value to variables
        pressure = self.get_p(coords[self.axial_dir]).reshape(-1, 1)
        interface_input.set_variable_data(self.mp_name_in, 'pressure', pressure)
        traction = np.zeros((model_part.size, 3))
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

        indices = sorted([self.axial_dir] + self.radial_dirs)  # columns that contain non-zero data
        a4_extra = np.delete(self.a4, indices, axis=1)  # remove columns containing data
        np.testing.assert_array_equal(a4_extra, a4_extra * 0.0)  # if a column remains it should be all zeroes
        np.testing.assert_allclose(self.a4[:, indices], self.a1[:, indices], rtol=1e-10, atol=1e-17)  # non-zero columns

    @unittest.skipIf(dimension == 2, 'Traction not supported in 2D')
    def test_shear(self):
        """
        Test whether shear is also applied.
        """

        # create solver
        solver = create_instance(self.parameters)
        solver.initialize()
        interface_input = solver.get_interface_input()
        model_part = interface_input.get_model_part(self.mp_name_in)
        coords = [model_part.x0, model_part.y0, model_part.z0]

        # give value to variables
        pressure = self.get_p(coords[self.axial_dir]).reshape(-1, 1)
        interface_input.set_variable_data(self.mp_name_in, 'pressure', pressure)
        traction = self.get_shear(coords[self.axial_dir], self.axial_dir)
        interface_input.set_variable_data(self.mp_name_in, 'traction', traction)

        # do 4 steps
        for i in range(4):
            solver.initialize_solution_step()
            solver.solve_solution_step(interface_input)
            solver.finalize_solution_step()
            solver.output_solution_step()
        solver.finalize()

        # compare output, as shear input has changed these should be different
        output_shear = solver.get_interface_output()
        self.a5 = output_shear.get_variable_data(self.mp_name_out, 'displacement')

        # calculate mean displacement
        # no absolute value is used on purpose to filter out the axial displacement due to pressure only
        self.mean_displacement_shear = np.mean(self.a5[:, self.axial_dir])
        self.mean_displacement_no_shear = np.mean(self.a1[:, self.axial_dir])

        print(f'Mean displacement in axial direction without shear = {self.mean_displacement_no_shear} m')
        print(f'Mean displacement in axial direction with shear = {self.mean_displacement_shear} m')
        self.assertNotAlmostEqual(self.mean_displacement_no_shear - self.mean_displacement_shear, 0., delta=1e-12)


class TestSolverWrapperAbaqusCSETube3D(TestSolverWrapperAbaqusCSETube2D):
    version = None
    setup_case = True
    dimension = 3
    axial_dir = 0  # x-direction is axial direction
    radial_dirs = [1, 2]
    mp_name_in = 'WALLOUTSIDE_load_points'
    mp_name_out = 'WALLOUTSIDE_nodes'


if __name__ == '__main__':
    unittest.main()
