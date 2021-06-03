from coconut.tools import create_instance, solver_available

import numpy as np
from os import getcwd
from os.path import realpath, dirname, join, exists
import unittest
import json
import pandas as pd
from shutil import copytree, rmtree


version = '60'


@unittest.skipUnless(solver_available(f'kratos.structural_mechanics_application.v{version}'),
                     f'kratos.structural_mechanics_application.v{version} not available')
class TestSolverWrapperKratosStructure60(unittest.TestCase):

    def setUp(self):
        folder_name = 'test_structure_tube_v60'
        dir_test = join(realpath(dirname(__file__)), folder_name)

        # Create setup files for kratos
        self.dir_csm = join(dir_test, 'CSM')
        if exists(self.dir_csm):
            rmtree(self.dir_csm)
        copytree(join(dir_test, 'setup_kratos'), join(dir_test, 'CSM'))

        parameter_file_name = join(dir_test, 'test_solver_wrapper.json')

        with open(parameter_file_name, 'r') as parameter_file:
            parameters = json.load(parameter_file)
            self.solver_param = parameters['solver_wrappers'][0]

        if getcwd() == realpath(dirname(__file__)):
            self.solver_param['settings']['working_directory'] = f'{folder_name}/CSM'

        self.mp_out_name = self.solver_param['settings']['interface_output'][0]['model_part']

    def tearDown(self):
        rmtree(self.dir_csm)

    # test if the same load results in the same displacement.
    def test_apply_load(self):
        solver = create_instance(self.solver_param)
        self.mp_out = solver.model.get_model_part(self.mp_out_name)

        load_interface = solver.get_interface_input()

        initial_pressure = 100

        pressure = initial_pressure
        traction_x = 0.0
        traction_y = 0.0
        traction_z = 0.0
        traction = np.array([traction_x, traction_y, traction_z])
        load_data_1 = self.get_uniform_load_data(pressure, traction)
        load_interface.set_interface_data(load_data_1)

        solver.initialize()

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
        solver.finalize()

        # obtain the  data and compare
        a1 = output_1.get_interface_data()
        a2 = output_2.get_interface_data()
        a3 = output_3.get_interface_data()

        load_data_to_kratos = self.read_pressure_traction()

        np.testing.assert_allclose(load_data_to_kratos, load_data_1, rtol=1e-12)
        np.testing.assert_allclose(a1, a3, rtol=0, atol=1e-10)
        self.assertRaises(AssertionError, np.testing.assert_allclose, a1, a2, rtol=0, atol=1e-10)

    # restart only works for memmabrane elements among planar elements (for kratos 6.0)
    def test_restart(self):
        self.solver_param['settings']['save_restart'] = -2
        solver = create_instance(self.solver_param)
        self.mp_out = solver.model.get_model_part(self.mp_out_name)
        load_interface = solver.get_interface_input()

        initial_pressure = 10

        pressure = initial_pressure
        traction_x = 0.0
        traction_y = 0.0
        traction_z = 0.0
        traction = np.array([traction_x, traction_y, traction_z])

        solver.initialize()

        # run solver for 4 timesteps
        for i in range(4):
            load = self.get_uniform_load_data(pressure * i, traction * i)
            load_interface.set_interface_data(load)
            solver.initialize_solution_step()
            solver.solve_solution_step(load_interface).copy()
            solver.finalize_solution_step()

        interface_x_1 = solver.get_interface_input()
        interface_y_1 = solver.get_interface_output()
        solver.finalize()

        # get data for solver without restart
        interface_output = solver.get_interface_output()
        out_data = interface_output.get_interface_data()

        self.solver_param['settings']['timestep_start'] = 2
        solver = create_instance(self.solver_param)
        solver.initialize()
        for i in range(2, 4):
            load = self.get_uniform_load_data(pressure * i, traction * i)
            load_interface.set_interface_data(load)
            solver.initialize_solution_step()
            solver.solve_solution_step(load_interface).copy()
            solver.finalize_solution_step()

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
        load = []
        for i in range(0, self.mp_out.size):
            load.append(pressure)
        for i in range(0, self.mp_out.size):
            load += traction.tolist()

        return np.array(load)

    def get_non_uniform_load(self, pressure, traction):
        l0 = 0.05
        pressure_data = pressure * np.sin(2 * np.pi / l0 * self.mp_out.x0)
        traction_x = traction[0] * np.sin(2 * np.pi / l0 * self.mp_out.x0)
        traction_y = traction[1] * np.sin(2 * np.pi / l0 * self.mp_out.x0)
        traction_z = traction[2] * np.sin(2 * np.pi / l0 * self.mp_out.x0)
        traction_data = np.column_stack((traction_x, traction_y, traction_z))
        traction_data = np.ravel(traction_data)

        return np.concatenate((pressure_data, traction_data))

    def read_pressure_traction(self):
        working_directory = self.solver_param['settings']['working_directory']
        input_mp_name = self.solver_param['settings']['kratos_interface_sub_model_parts_list'][0]
        pr_fname = join(working_directory, f'{input_mp_name}_pressure.csv')
        sl_fname = join(working_directory, f'{input_mp_name}_surface_load.csv')
        df_pr = pd.read_csv(pr_fname, skipinitialspace=True)
        df_sl = pd.read_csv(sl_fname, skipinitialspace=True)
        pressure = np.array(df_pr['pressure'])
        traction_x = np.array(df_sl['surface_load_x'])
        traction_y = np.array(df_sl['surface_load_y'])
        traction_z = np.array(df_sl['surface_load_z'])
        traction = np.ravel(np.column_stack((traction_x, traction_y, traction_z)))
        return np.concatenate((pressure, traction))


if __name__ == '__main__':
    unittest.main()
