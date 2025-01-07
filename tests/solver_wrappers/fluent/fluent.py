from coconut.tools import create_instance, get_solver_env, rm_timed

import unittest
import numpy as np
import os
from os.path import join
import subprocess
import json
import multiprocessing
import shutil


class TestSolverWrapperFluentTube2D(unittest.TestCase):
    version = 'xxxxRx'  # Fluent product version, as from 2023R1 typically of the form 'xxxRx', set in subclass
    setup_case = True

    @classmethod
    def setUpClass(cls):
        dimension = 2
        dir_name = os.path.realpath(os.path.dirname(__file__))  # path to fluent directory
        cls.file_name = join(dir_name, f'test_v{cls.version}/tube{dimension}d/parameters.json')
        cls.working_dir = join(dir_name, f'test_v{cls.version}/tube{dimension}d/CFD')
        cls.setup_dir = join(dir_name, f'test_v{cls.version}/tube{dimension}d/setup_fluent')

        # setup
        if cls.setup_case:
            shutil.rmtree(cls.working_dir, ignore_errors=True)
            shutil.copytree(cls.setup_dir, cls.working_dir)
            fluent_solver_module = f'fluent.v{cls.version}'
            env = get_solver_env(fluent_solver_module, cls.working_dir)
            subprocess.check_call(join(cls.working_dir, 'setup_fluent.sh'), shell=True, cwd=cls.working_dir, env=env)

    def setUp(self):
        with open(self.file_name) as parameter_file:
            self.parameters = json.load(parameter_file)
        self.parameters['settings']['working_directory'] = os.path.relpath(self.working_dir)  # set working directory
        self.parameters['settings']['cores'] = min(4, multiprocessing.cpu_count())
        self.mp_name_in = 'beamoutside_nodes'
        self.mp_name_out = 'beamoutside_faces'

    #@classmethod
    #def tearDownClass(cls):
    #    if cls.setup_case:
    #        rm_timed(cls.working_dir)

    # noinspection PyMethodMayBeStatic
    def get_dy(self, x):
        return 0.0001 * np.sin(2 * np.pi / 0.05 * x)

    def test_move_nodes(self):
        # test if nodes are moved to the correct position

        # adapt parameters, create solver
        self.parameters['settings']['flow_iterations'] = 1
        solver = create_instance(self.parameters)
        solver.initialize()

        # set displacement
        interface_input = solver.get_interface_input()
        model_part = interface_input.get_model_part(self.mp_name_in)

        x0, y0 = model_part.x0, model_part.y0
        dy = self.get_dy(x0)

        displacement = interface_input.get_variable_data(self.mp_name_in, 'displacement')
        displacement[:, 1] = dy
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
        np.testing.assert_allclose(y, y0 + dy, rtol=1e-15)

    def test_partitioning(self):
        # test if different partitioning gives the same ModelParts

        # create two solvers with different partitioning
        x0, y0, z0, ids = [], [], [], []
        for cores in [1, max(2, min(4, multiprocessing.cpu_count()))]:
            self.parameters['settings']['cores'] = cores
            solver = create_instance(self.parameters)
            solver.initialize()
            solver.finalize()
            mp = solver.model.get_model_part(self.mp_name_in)
            x0.append(mp.x0)
            y0.append(mp.y0)
            z0.append(mp.z0)
            ids.append(mp.id)

        # compare ModelParts of both solvers
        for attr in [x0, y0, z0, ids]:
            np.testing.assert_array_equal(attr[0], attr[1])

    def test_pressure_traction(self):
        # test if same coordinates always give same pressure & traction

        self.parameters['settings']['cores'] = min(4, multiprocessing.cpu_count())
        self.parameters['settings']['flow_iterations'] = 500
        solver = create_instance(self.parameters)
        solver.initialize()
        solver.initialize_solution_step()
        interface_input = solver.get_interface_input()

        # set displacement
        x0 = interface_input.get_model_part(self.mp_name_in).x0
        dys = [self.get_dy(x0), np.zeros_like(x0), self.get_dy(x0)]
        displacement = interface_input.get_variable_data(self.mp_name_in, 'displacement')

        # run solver for three displacements (first one = last one)
        pressure = []
        traction = []
        for dy in dys:
            displacement[:, 1] = dy
            interface_input.set_variable_data(self.mp_name_in, 'displacement', displacement)
            interface_output = solver.solve_solution_step(interface_input)
            pressure.append(interface_output.get_variable_data(self.mp_name_out, 'pressure'))
            traction.append(interface_output.get_variable_data(self.mp_name_out, 'traction'))
        solver.finalize_solution_step()
        solver.output_solution_step()
        solver.finalize()

        # check if same position gives same pressure & traction
        np.testing.assert_allclose(pressure[0], pressure[2], rtol=1e-10)
        np.testing.assert_allclose(traction[0], traction[2], rtol=1e-11)

        # check if different position gives different pressure & traction
        p01 = np.linalg.norm(pressure[0] - pressure[1])
        p02 = np.linalg.norm(pressure[0] - pressure[2])
        self.assertTrue(p02 / p01 < 1e-12)

        t01 = np.linalg.norm(traction[0] - traction[1])
        t02 = np.linalg.norm(traction[0] - traction[2])
        self.assertTrue(t02 / t01 < 1e-12)

    def test_restart(self):
        # test if restart option works correctly

        # adapt parameters, create solver
        self.parameters['settings']['cores'] = 1
        self.parameters['settings']['flow_iterations'] = 30
        solver = create_instance(self.parameters)
        solver.initialize()
        interface_input = solver.get_interface_input()

        # set displacement
        x0 = interface_input.get_model_part(self.mp_name_in).x0
        displacement = interface_input.get_variable_data(self.mp_name_in, 'displacement')
        displacement[:, 1] = self.get_dy(x0)

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
        pressure_1 = interface_output.get_variable_data(self.mp_name_out, 'pressure')
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
        pressure_2 = interface_output.get_variable_data(self.mp_name_out, 'pressure')
        traction_2 = interface_output.get_variable_data(self.mp_name_out, 'traction')

        # check if undeformed coordinate (coordinates of model part) are equal
        self.assertTrue(interface_x_1.has_same_model_parts(interface_x_2))
        self.assertTrue(interface_y_1.has_same_model_parts(interface_y_2))

        # check if coordinates of ModelParts are equal
        # ==>  check if deformed coordinates are equal
        np.testing.assert_allclose(coords_1, coords_2, rtol=1e-15)

        # check if pressure and traction are equal
        np.testing.assert_allclose(pressure_1, pressure_2, rtol=1e-14)
        np.testing.assert_allclose(traction_1, traction_2, rtol=1e-14)

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
        x0 = interface_input.get_model_part(self.mp_name_in).x0
        displacement = interface_input.get_variable_data(self.mp_name_in, 'displacement')
        displacement[:, 1] = self.get_dy(x0)

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


class TestSolverWrapperFluentTube3D(unittest.TestCase):
    version = 'xxxxRx'  # Fluent product version, as from 2023R1 typically of the form 'xxxRx', set in subclass
    setup_case = True

    @classmethod
    def setUpClass(cls):
        dimension = 3
        dir_name = os.path.realpath(os.path.dirname(__file__))  # path to fluent directory
        cls.file_name = join(dir_name, f'test_v{cls.version}/tube{dimension}d/parameters.json')
        cls.working_dir = join(dir_name, f'test_v{cls.version}/tube{dimension}d/CFD')
        cls.setup_dir = join(dir_name, f'test_v{cls.version}/tube{dimension}d/setup_fluent')

        # setup
        if cls.setup_case:
            shutil.rmtree(cls.working_dir, ignore_errors=True)
            shutil.copytree(cls.setup_dir, cls.working_dir)
            fluent_solver_module = f'fluent.v{cls.version}'
            env = get_solver_env(fluent_solver_module, cls.working_dir)
            subprocess.check_call(join(cls.working_dir, 'setup_fluent.sh'), shell=True, cwd=cls.working_dir, env=env)

    def setUp(self):
        with open(self.file_name) as parameter_file:
            self.parameters = json.load(parameter_file)
        self.parameters['settings']['working_directory'] = os.path.relpath(self.working_dir)  # set working directory
        self.parameters['settings']['cores'] = min(4, multiprocessing.cpu_count())
        self.mp_name_in = 'wall_nodes'
        self.mp_name_out = 'wall_faces'

    @classmethod
    def tearDownClass(cls):
        if cls.setup_case:
            rm_timed(cls.working_dir)

    # noinspection PyMethodMayBeStatic
    def get_dy_dz(self, x, y, z):
        dr = 0.0001 * np.sin(2 * np.pi / 0.05 * x)
        theta = np.arctan2(z, y)
        dy = dr * np.cos(theta)
        dz = dr * np.sin(theta)
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
        dy, dz = self.get_dy_dz(x0, y0, z0)

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
        # test if different partitioning gives the same ModelParts

        # create two solvers with different partitioning
        x0, y0, z0, ids = [], [], [], []
        for cores in [1, max(2, min(4, multiprocessing.cpu_count()))]:
            self.parameters['settings']['cores'] = cores
            solver = create_instance(self.parameters)
            solver.initialize()
            solver.finalize()
            mp = solver.model.get_model_part('wall_nodes')
            x0.append(mp.x0)
            y0.append(mp.y0)
            z0.append(mp.z0)
            ids.append(mp.id)

        # compare ModelParts of both solvers
        for attr in [x0, y0, z0, ids]:
            np.testing.assert_array_equal(attr[0], attr[1])

    def test_pressure_traction(self):
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
        dy, dz = self.get_dy_dz(x0, y0, z0)
        dydz = np.column_stack((dy, dz))
        dydzs = [dydz, np.zeros_like(dydz), dydz]
        displacement = interface_input.get_variable_data(self.mp_name_in, 'displacement')

        # run solver for three displacements (first one = last one)
        pressure = []
        traction = []
        for dydz in dydzs:
            displacement[:, 1:] = dydz
            interface_input.set_variable_data(self.mp_name_in, 'displacement', displacement)
            interface_output = solver.solve_solution_step(interface_input)
            pressure.append(interface_output.get_variable_data(self.mp_name_out, 'pressure'))
            traction.append(interface_output.get_variable_data(self.mp_name_out, 'traction'))
        solver.finalize_solution_step()
        solver.output_solution_step()
        solver.finalize()

        # check if same position gives same pressure & traction
        np.testing.assert_allclose(pressure[0], pressure[2], rtol=1e-10)
        np.testing.assert_allclose(traction[0], traction[2], rtol=1e-10)

        # check if different position gives different pressure & traction
        p01 = np.linalg.norm(pressure[0] - pressure[1])
        p02 = np.linalg.norm(pressure[0] - pressure[2])
        self.assertTrue(p02 / p01 < 1e-12)

        t01 = np.linalg.norm(traction[0] - traction[1])
        t02 = np.linalg.norm(traction[0] - traction[2])
        self.assertTrue(t02 / t01 < 1e-12)

    def test_restart(self):
        # test if restart option works correctly

        # adapt parameters, create solver
        self.parameters['settings']['cores'] = 1
        self.parameters['settings']['flow_iterations'] = 30
        solver = create_instance(self.parameters)
        solver.initialize()
        interface_input = solver.get_interface_input()

        # set displacement
        model_part = interface_input.get_model_part(self.mp_name_in)
        x0, y0, z0 = model_part.x0, model_part.y0, model_part.z0
        dy, dz = self.get_dy_dz(x0, y0, z0)
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
        pressure_1 = interface_output.get_variable_data(self.mp_name_out, 'pressure')
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
        pressure_2 = interface_output.get_variable_data(self.mp_name_out, 'pressure')
        traction_2 = interface_output.get_variable_data(self.mp_name_out, 'traction')

        # check if undeformed coordinate (coordinates of model part) are equal
        self.assertTrue(interface_x_1.has_same_model_parts(interface_x_2))
        self.assertTrue(interface_y_1.has_same_model_parts(interface_y_2))

        # check if coordinates of ModelParts are equal
        # ==>  check if deformed coordinates are equal
        np.testing.assert_allclose(coords_1, coords_2, rtol=1e-15)

        # check if pressure and traction are equal
        np.testing.assert_allclose(pressure_1, pressure_2, rtol=1e-14)
        np.testing.assert_allclose(traction_1, traction_2, rtol=1e-14)

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
        dy, dz = self.get_dy_dz(x0, y0, z0)
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


if __name__ == '__main__':
    unittest.main()
