from coconut import data_structure
from coconut.tools import create_instance
from coconut.data_structure.interface import Interface

import numpy as np
import os
import subprocess
import unittest
import json
import pandas as pd


def print_box(text):
    n = len(text)
    top = '\n┌─' + n * '─' + '─┐'
    mid = '\n│ ' + text + ' │'
    bottom = '\n└─' + n * '─' + '─┘'
    print(top + mid + bottom)


class TestSolverWrapperKratosStructure60(unittest.TestCase):

    def setUp(self):
        folder_name = 'test_structure_tube_v60'
        self.dir_tmp = os.path.join(os.path.realpath(os.path.dirname(__file__)), folder_name)
        setup_file = os.path.join(self.dir_tmp, 'setup_kratos.sh')
        p = subprocess.Popen(f'sh {setup_file}', cwd=self.dir_tmp, shell=True)
        p.wait()

        # Create setup files for kratos

        parameter_file_name = os.path.join(self.dir_tmp, 'test_solver_wrapper.json')

        with open(parameter_file_name, 'r') as parameter_file:
            parameters = json.load(parameter_file)
            self.solver_param = parameters['solver_wrappers'][0]

        if os.getcwd() == os.path.realpath(os.path.dirname(__file__)):
            self.solver_param['settings']['working_directory'] = f'{folder_name}/CSM'

        self.solver = create_instance(self.solver_param)
        mp_out_name = self.solver_param['settings']['interface_output'][0]['model_part']
        self.mp_out = self.solver.model.get_model_part(mp_out_name)

    def test_apply_load(self):
        print_box('started tests for Kratos Tube3D: test_apply_load')

        load_interface = self.solver.get_interface_input()

        initial_pressure = 100

        pressure = initial_pressure
        traction_x = 100.0
        traction_y = 0.0
        traction_z = 0.0
        traction = np.array([traction_x, traction_y, traction_z])
        load_data_1 = self.get_non_uniform_load(pressure, traction)
        load_interface.set_interface_data(load_data_1)

        self.solver.initialize()

        self.solver.initialize_solution_step()
        output_1 = self.solver.solve_solution_step(load_interface).copy()

        pressure *= 0.5
        traction *= 0.5
        load_data_2 = self.get_non_uniform_load(pressure, traction)
        load_interface.set_interface_data(load_data_2)
        output_2 = self.solver.solve_solution_step(load_interface).copy()

        pressure = initial_pressure
        traction = np.array([traction_x, traction_y, traction_z])
        load_data_3 = self.get_non_uniform_load(pressure, traction)
        load_interface.set_interface_data(load_data_3)
        output_3 = self.solver.solve_solution_step(load_interface).copy()
        self.solver.finalize_solution_step()
        self.solver.finalize()

        # obtain the  data and compare
        a1 = output_1.get_interface_data()
        a2 = output_2.get_interface_data()
        a3 = output_3.get_interface_data()

        load_data_to_kratos = self.read_pressure_traction()

        np.testing.assert_allclose(load_data_to_kratos, load_data_1, rtol=1e-12)
        np.testing.assert_allclose(a1, a3, rtol=0, atol=1e-12)
        self.assertRaises(AssertionError, np.testing.assert_allclose, a1, a2, rtol=0, atol=1e-12)

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
        pr_fname = os.path.join(working_directory, f'{input_mp_name}_pressure.csv')
        sl_fname = os.path.join(working_directory, f'{input_mp_name}_surface_load.csv')
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
