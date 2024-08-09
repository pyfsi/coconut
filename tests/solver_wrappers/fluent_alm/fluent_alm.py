from coconut.tools import create_instance, get_solver_env, rm_timed

import unittest
import numpy as np
import os
from os.path import join
import subprocess
import json
import multiprocessing
import shutil
from scipy import integrate


class TestSolverWrapperFluentALMYarn(unittest.TestCase):
    version = 'xxxxRx'  # Fluent product version, as from 2023R1 typically of the form 'xxxRx', set in subclass
    setup_case = True

    @classmethod
    def setUpClass(cls):
        dir_name = os.path.realpath(os.path.dirname(__file__))  # path to alm_fluent directory
        cls.file_name = join(dir_name, f'test_v{cls.version}/yarn3d/parameters.json')
        cls.working_dir = join(dir_name, f'test_v{cls.version}/yarn3d/CFD')
        # setup
        if cls.setup_case:
            shutil.rmtree(cls.working_dir, ignore_errors=True)
            dir_tmp = join(dir_name, f'test_v{cls.version}/yarn3d')
            fluent_solver_module = f'fluent_alm.v{cls.version}'
            env = get_solver_env(fluent_solver_module, dir_tmp)
            p = subprocess.Popen(join(dir_tmp, 'setup_fluent.sh'), cwd=dir_tmp, shell=True, env=env)
            p.wait()

    def setUp(self):
        with open(self.file_name) as parameter_file:
            self.parameters = json.load(parameter_file)
        self.parameters['settings']['working_directory'] = os.path.relpath(self.working_dir)  # set working directory
        self.parameters['settings']['cores'] = min(4, multiprocessing.cpu_count())
        self.mp_name_in = 'yarn_coordinate_points'
        self.mp_name_out = 'yarn_load_points'

    @classmethod
    def tearDownClass(cls):
        if cls.setup_case:
            rm_timed(cls.working_dir)

    # noinspection PyMethodMayBeStatic
    def get_dy_dz(self, x):
        r = np.abs(x)
        theta = np.radians(1.)
        dy = 0.005 * np.sin((x - x[0]) / 0.25 * np.pi)
        dz = r * np.sin(theta)
        return dy, dz

    def test_move_nodes(self):
        # test if nodes are moved to the correct position

        # adapt parameters, create solver
        self.parameters['settings']['flow_iterations'] = 1
        solver = create_instance(self.parameters)
        solver.initialize()

        # set displacement
        interface_input = solver.get_interface_input()
        model_part = interface_input.get_model_part(self.mp_name_in)

        x0, y0, z0 = model_part.x0, model_part.y0, model_part.z0
        dy, dz = self.get_dy_dz(x0)

        displacement = interface_input.get_variable_data(self.mp_name_in, 'displacement')
        displacement[:, 1] = dy
        displacement[:, 2] = dz
        interface_input.set_variable_data(self.mp_name_in, 'displacement', displacement)

        # update position by iterating once in solver
        solver.initialize_solution_step()
        solver.solve_solution_step(interface_input)
        solver.finalize_solution_step()
        solver.output_solution_step()
        coord_data = solver.get_coordinates()
        solver.finalize()

        # check if correct displacement was given
        y = coord_data[self.mp_name_in]['coords'][:, 1]
        z = coord_data[self.mp_name_in]['coords'][:, 2]
        np.testing.assert_allclose(y, y0 + dy, rtol=1e-15)
        np.testing.assert_allclose(z, z0 + dz, rtol=1e-15)

    def test_partitioning(self):
        # test if different partitioning gives the same aerodynamic forces

        """
        NOTE: for the partitioning test, it is important that the actuator points do NOT coincide with mesh points.
        If this happens, the find_cell_at_yarn_points-function in the UDF can/will assign different cells as 'owner'
        for this yarn points, for different partitioning. This then causes a slight change in sampled density (as this
        is a cell-centered sampling), leading to a change in aerodynamic forces and eventually velocities as well.
        I noticed for this unit test that the relative difference is then in the order of 1e-4, causing the test to fail.
        """

        self.parameters['settings']['flow_iterations'] = 250

        # create two solvers with different partitioning and solve solution step
        traction = []
        for cores in [1, max(2, min(4, multiprocessing.cpu_count()))]:
            # create solver
            self.parameters['settings']['cores'] = cores
            solver = create_instance(self.parameters)
            solver.initialize()
            solver.initialize_solution_step()
            interface_input = solver.get_interface_input()

            # set displacement
            model_part = interface_input.get_model_part(self.mp_name_in)
            x0, y0, z0 = model_part.x0, model_part.y0, model_part.z0
            dy, dz = self.get_dy_dz(x0)
            displacement = interface_input.get_variable_data(self.mp_name_in, 'displacement')
            displacement[:, 1] = dy
            displacement[:, 2] = dz
            interface_input.set_variable_data(self.mp_name_in, 'displacement', displacement)

            # solve solution step and output data
            interface_output = solver.solve_solution_step(interface_input)
            traction.append(interface_output.get_variable_data(self.mp_name_out, 'traction'))

            # finalize solver
            solver.finalize_solution_step()
            solver.output_solution_step()
            solver.finalize()

        # compare output traction of both solvers
        np.testing.assert_allclose(traction[0], traction[1], rtol=1e-10, atol=1e-12)

    def test_traction(self):
        # test if same coordinates always give same pressure & traction

        self.parameters['settings']['cores'] = min(4, multiprocessing.cpu_count())
        self.parameters['settings']['flow_iterations'] = 500
        solver = create_instance(self.parameters)
        solver.initialize()
        solver.initialize_solution_step()
        interface_input = solver.get_interface_input()

        # set displacement
        model_part = interface_input.get_model_part(self.mp_name_in)
        x0, y0, z0 = model_part.x0, model_part.y0, model_part.z0
        dy, dz = self.get_dy_dz(x0)
        dydz = np.column_stack((dy, dz))
        dydzs = [dydz, np.zeros_like(dydz), dydz]
        displacement = interface_input.get_variable_data(self.mp_name_in, 'displacement')

        # run solver for three displacements (first one = last one)
        traction = []
        for dydz in dydzs:
            displacement[:, 1:] = dydz
            interface_input.set_variable_data(self.mp_name_in, 'displacement', displacement)
            interface_output = solver.solve_solution_step(interface_input)
            traction.append(interface_output.get_variable_data(self.mp_name_out, 'traction'))
        solver.finalize_solution_step()
        solver.output_solution_step()
        solver.finalize()

        # check if same position gives same traction
        np.testing.assert_allclose(traction[0], traction[2], rtol=1e-10)

        # check if different position gives different traction
        t01 = np.linalg.norm(traction[0] - traction[1])
        t02 = np.linalg.norm(traction[0] - traction[2])
        self.assertTrue(t02 / t01 < 1e-12)

    def test_restart(self):
        # test if restart option works correctly

        # adapt parameters, create solver
        self.parameters['settings']['cores'] = 1
        self.parameters['settings']['flow_iterations'] = 50
        solver = create_instance(self.parameters)
        solver.initialize()
        interface_input = solver.get_interface_input()

        # set displacement
        model_part = interface_input.get_model_part(self.mp_name_in)
        x0, y0, z0 = model_part.x0, model_part.y0, model_part.z0
        dy, dz = self.get_dy_dz(x0)
        displacement = interface_input.get_variable_data(self.mp_name_in, 'displacement')
        displacement[:, 1] = dy
        displacement[:, 2] = dz

        # run solver for 4 time steps
        for i in range(4):
            solver.initialize_solution_step()
            interface_input.set_variable_data(self.mp_name_in, 'displacement', i * displacement)
            solver.solve_solution_step(interface_input)
            solver.finalize_solution_step()
            solver.output_solution_step()
        interface_x_1 = solver.get_interface_input()
        interface_y_1 = solver.get_interface_output()
        coords_1 = solver.get_coordinates()[self.mp_name_in]['coords']
        solver.finalize()

        # get data for solver without restart
        interface_output = solver.get_interface_output()
        traction_1 = interface_output.get_variable_data(self.mp_name_out, 'traction')

        # create solver which restarts at time step 2
        self.parameters['settings']['timestep_start'] = 2
        solver = create_instance(self.parameters)
        solver.initialize()
        interface_input = solver.get_interface_input()

        # run solver for 2 more time steps
        for i in range(2, 4):
            solver.initialize_solution_step()
            interface_input.set_variable_data(self.mp_name_in, 'displacement', i * displacement)
            solver.solve_solution_step(interface_input)
            solver.finalize_solution_step()
            solver.output_solution_step()
        interface_x_2 = solver.get_interface_input()
        interface_y_2 = solver.get_interface_output()
        coords_2 = solver.get_coordinates()[self.mp_name_in]['coords']
        solver.finalize()

        # get data for solver with restart
        interface_output = solver.get_interface_output()
        traction_2 = interface_output.get_variable_data(self.mp_name_out, 'traction')

        # check if undeformed coordinate (coordinates of model part) are equal
        self.assertTrue(interface_x_1.has_same_model_parts(interface_x_2))
        self.assertTrue(interface_y_1.has_same_model_parts(interface_y_2))

        # check if coordinates of ModelParts are equal
        # ==>  check if deformed coordinates are equal
        np.testing.assert_allclose(coords_1, coords_2, rtol=1e-15)

        # check if traction is equal
        np.testing.assert_allclose(traction_1, traction_2, rtol=1e-12)

    def test_coupling_convergence(self):
        # test if check of coupling convergence works correctly

        # adapt parameters, create solver
        self.parameters['settings']['cores'] = 1
        self.parameters['settings']['flow_iterations'] = 100
        solver = create_instance(self.parameters)
        solver.check_coupling_convergence = True
        solver.initialize()
        interface_input = solver.get_interface_input()

        # set displacement
        model_part = interface_input.get_model_part(self.mp_name_in)
        x0, y0, z0 = model_part.x0, model_part.y0, model_part.z0
        dy, dz = self.get_dy_dz(x0)
        displacement = interface_input.get_variable_data(self.mp_name_in, 'displacement')
        displacement[:, 1] = dy
        displacement[:, 2] = dz

        solver.initialize_solution_step()

        # first coupling iteration
        interface_input.set_variable_data(self.mp_name_in, 'displacement', displacement)
        solver.solve_solution_step(interface_input)

        self.assertFalse(solver.coupling_convergence)

        # second coupling iteration
        interface_input.set_variable_data(self.mp_name_in, 'displacement', 2 * displacement)
        solver.solve_solution_step(interface_input)

        self.assertFalse(solver.coupling_convergence)

        # third coupling iteration
        interface_input.set_variable_data(self.mp_name_in, 'displacement', 2 * displacement)
        solver.solve_solution_step(interface_input)

        self.assertTrue(solver.coupling_convergence)

        solver.output_solution_step()
        solver.finalize_solution_step()
        solver.finalize()

    def test_forces(self):
        # test if volume integral of momentum sources is nearly equal to line integral of actuator forces

        # remove all previous 'report-coconut*'-files so that 'report-coconut.out' is free for Fluent to write
        for fname in os.listdir(self.working_dir):
            if fname.startswith('report-coconut'):
                os.remove(join(self.working_dir, fname))

        # create solver
        self.parameters['settings']['flow_iterations'] = 50
        solver = create_instance(self.parameters)
        solver.initialize()
        interface_input = solver.get_interface_input()
        solver.initialize_solution_step()

        # set zero displacement (initial position: aligned with flow, constant actuator point spacing)
        displacement = interface_input.get_variable_data(self.mp_name_in, 'displacement') * 0.
        interface_input.set_variable_data(self.mp_name_in, 'displacement', displacement)

        # solve solution step and output data
        interface_output = solver.solve_solution_step(interface_input)
        traction = interface_output.get_variable_data(self.mp_name_out, 'traction')
        line_force = integrate.trapezoid(traction,
                                         dx=0.5 / (interface_output.get_model_part(self.mp_name_out).size - 1), axis=0)

        # finalize solver
        solver.finalize_solution_step()
        solver.output_solution_step()
        solver.finalize()

        # get volumetric forces from report file
        report_file = open(join(self.working_dir, 'report-coconut.out'), 'r')
        last_line = report_file.readlines()[-1]
        volume_force = -1 * np.array(last_line.strip().split()[1:4], dtype=float)
        report_file.close()

        # test if force vectors are almost equal (up to 1.5 %)
        np.testing.assert_allclose(line_force, volume_force, rtol=1.5e-2)


if __name__ == '__main__':
    unittest.main()
