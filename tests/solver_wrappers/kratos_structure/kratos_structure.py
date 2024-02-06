from coconut.tools import create_instance

import numpy as np
from os.path import realpath, dirname, join, relpath
import unittest
import json
import pandas as pd
from shutil import copytree, rmtree


class BaseTestSolverWrapperKratosStructure(unittest.TestCase):
    version_label = None  # KratosMultiphysics version without dot, e.g. 91, set in subclass

    @classmethod
    def setUpClass(cls):
        dir_name = join(realpath(dirname(__file__)))
        dir_test = f'test_v{cls.version_label}'

        # create setup files for kratos
        cls.working_dir = join(dir_name, dir_test, 'CSM')
        rmtree(cls.working_dir, ignore_errors=True)
        copytree(join(dir_name, dir_test, 'setup_kratos'), cls.working_dir)
        cls.parameter_file_name = join(dir_name, dir_test, 'parameters.json')

    def setUp(self):
        with open(self.parameter_file_name, 'r') as parameter_file:
            self.parameters = json.load(parameter_file)

        self.parameters['settings']['working_directory'] = relpath(self.working_dir)

        self.mp_in_name = self.parameters['settings']['interface_input'][0]['model_part']

    @classmethod
    def tearDownClass(cls):
        rmtree(cls.working_dir, ignore_errors=True)

    # test if the same load results in the same displacement.
    def test_apply_load(self):
        solver = create_instance(self.parameters)
        solver.initialize()
        self.mp_in = solver.model.get_model_part(self.mp_in_name)

        load_interface = solver.get_interface_input()

        initial_pressure = 100

        pressure = initial_pressure
        traction_x = 0.3  # set traction values to zero for faster testing
        traction_y = 0.2
        traction_z = 0.1
        traction = np.array([traction_x, traction_y, traction_z])
        load_data_1 = self.get_uniform_load_data(pressure, traction)
        load_interface.set_interface_data(load_data_1)

        solver.initialize_solution_step()
        output_1 = solver.solve_solution_step(load_interface).copy()

        pressure *= 0.5
        traction *= 0.5
        load_data_2 = self.get_uniform_load_data(pressure, traction)
        load_interface.set_interface_data(load_data_2)
        output_2 = solver.solve_solution_step(load_interface).copy()

        pressure = initial_pressure
        traction = np.array([traction_x, traction_y, traction_z])
        load_data_3 = self.get_uniform_load_data(pressure, traction)
        load_interface.set_interface_data(load_data_3)
        output_3 = solver.solve_solution_step(load_interface).copy()
        solver.finalize_solution_step()
        solver.output_solution_step()
        solver.finalize()

        # obtain the  data and compare
        a1 = output_1.get_interface_data()
        a2 = output_2.get_interface_data()
        a3 = output_3.get_interface_data()

        load_data_to_kratos = self.read_pressure_traction()

        np.testing.assert_allclose(load_data_to_kratos, load_data_1, rtol=1e-12)
        np.testing.assert_allclose(a1, a3, rtol=0, atol=1e-10)
        self.assertRaises(AssertionError, np.testing.assert_allclose, a1, a2, rtol=0, atol=1e-10)

    def test_restart(self):
        self.parameters['settings']['save_restart'] = 2
        solver = create_instance(self.parameters)
        solver.initialize()
        self.mp_in = solver.model.get_model_part(self.mp_in_name)
        load_interface = solver.get_interface_input()

        initial_pressure = 10

        pressure = initial_pressure
        traction_x = 0.0
        traction_y = 0.0
        traction_z = 0.0
        traction = np.array([traction_x, traction_y, traction_z])

        # run solver for 4 time steps
        for i in range(4):
            load = self.get_uniform_load_data(pressure * i, traction * i)
            load_interface.set_interface_data(load)
            solver.initialize_solution_step()
            solver.solve_solution_step(load_interface).copy()
            solver.finalize_solution_step()
            solver.output_solution_step()

        interface_x_1 = solver.get_interface_input()
        interface_y_1 = solver.get_interface_output()
        solver.finalize()

        # get data for solver without restart
        interface_output = solver.get_interface_output()
        out_data = interface_output.get_interface_data()

        self.parameters['settings']['timestep_start'] = 2
        solver = create_instance(self.parameters)
        solver.initialize()
        for i in range(2, 4):
            load = self.get_uniform_load_data(pressure * i, traction * i)
            load_interface.set_interface_data(load)
            solver.initialize_solution_step()
            solver.solve_solution_step(load_interface).copy()
            solver.finalize_solution_step()
            solver.output_solution_step()

        interface_x_2 = solver.get_interface_input()
        interface_y_2 = solver.get_interface_output()

        solver.finalize()

        # get data for solver with restart
        interface_output_restart = solver.get_interface_output()
        out_data_restart = interface_output_restart.get_interface_data()

        # check if undeformed coordinate (coordinates of model part) are equal
        self.assertTrue(interface_x_1.has_same_model_parts(interface_x_2))
        self.assertTrue(interface_y_1.has_same_model_parts(interface_y_2))

        max_value = np.max(np.abs(out_data))

        # check if pressure and traction are equal
        np.testing.assert_allclose(out_data / max_value, out_data_restart / max_value, rtol=0, atol=1e-10)

    def get_uniform_load_data(self, pressure, traction):
        load = self.mp_in.size * [pressure] + self.mp_in.size * traction.tolist()
        return np.array(load)

    def get_non_uniform_load_data(self, pressure, traction):
        l0 = 0.05
        pressure_data = pressure * np.sin(2 * np.pi / l0 * self.mp_in.x0)
        traction_x = traction[0] * np.sin(2 * np.pi / l0 * self.mp_in.x0)
        traction_y = traction[1] * np.sin(2 * np.pi / l0 * self.mp_in.x0)
        traction_z = traction[2] * np.sin(2 * np.pi / l0 * self.mp_in.x0)
        traction_data = np.column_stack((traction_x, traction_y, traction_z))
        traction_data = np.ravel(traction_data)

        return np.concatenate((pressure_data, traction_data))

    def read_pressure_traction(self):
        working_directory = self.parameters['settings']['working_directory']
        input_mp_name = self.parameters['settings']['kratos_interface_sub_model_parts_list'][0]
        pr_fname = join(working_directory, f'{input_mp_name}_pressure.csv')
        tr_fname = join(working_directory, f'{input_mp_name}_traction.csv')
        df_pr = pd.read_csv(pr_fname, skipinitialspace=True)
        df_tr = pd.read_csv(tr_fname, skipinitialspace=True)
        pressure = np.array(df_pr['pressure'])
        traction_x = np.array(df_tr['traction_x'])
        traction_y = np.array(df_tr['traction_y'])
        traction_z = np.array(df_tr['traction_z'])
        traction = np.ravel(np.column_stack((traction_x, traction_y, traction_z)))
        return np.concatenate((pressure, traction))

    def test_coupling_convergence(self):
        # test if check of coupling convergence works correctly

        # adapt parameters, create solver
        solver = create_instance(self.parameters)
        solver.check_coupling_convergence = True
        solver.initialize()
        self.mp_in = solver.model.get_model_part(self.mp_in_name)
        interface_input = solver.get_interface_input()

        # set pressure and traction
        initial_pressure = 100

        pressure = initial_pressure
        traction_x = 0.3  # set traction values to zero for faster testing
        traction_y = 0.2
        traction_z = 0.1
        traction = np.array([traction_x, traction_y, traction_z])
        load_data = self.get_non_uniform_load_data(pressure, traction)

        solver.initialize_solution_step()

        # first coupling iteration
        interface_input.set_interface_data(load_data)
        solver.solve_solution_step(interface_input)

        self.assertFalse(solver.coupling_convergence)

        # second coupling iteration
        interface_input.set_interface_data(load_data * 2)
        solver.solve_solution_step(interface_input)

        self.assertFalse(solver.coupling_convergence)

        # third coupling iteration
        interface_input.set_interface_data(load_data * 2)
        solver.solve_solution_step(interface_input)

        self.assertTrue(solver.coupling_convergence)

        solver.output_solution_step()
        solver.finalize_solution_step()
        solver.finalize()
