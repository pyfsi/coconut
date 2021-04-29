from coconut.tools import create_instance

import unittest
import os
import json
import numpy as np
from unittest import mock

def mock_create_instance(settings):
    object_type = settings['type']
    if not 'dummy' in object_type:
        return create_instance(settings)
    else:
        object_module = __import__('coconut.tests.solver_wrappers.combined.' + object_type, fromlist=[object_type])
    return object_module.create(settings)


class TestSolverWrapperCombined(unittest.TestCase):

    def setUp(self):

        parameter_file_name = os.path.join(os.path.dirname(__file__), 'parameters.json')

        with open(parameter_file_name, 'r') as parameter_file:
            parameters = json.load(parameter_file)
        self.master_index = parameters['settings']['master_solver_index']
        self.master_par_solver = parameters['settings']['solver_wrappers'][self.master_index]

        self.comb_sol_par = parameters

        self.mp_name_in = self.master_par_solver['settings']['interface_input'][0]['model_part']
        self.mp_name_out = self.master_par_solver['settings']['interface_output'][0]['model_part']


    def get_displacement(self, x, y, z):
        displacement = np.empty((x.size, 3))
        displacement[:,0] = 0.0
        displacement[:, 1] = 0.0
        displacement[:, 2] = (x[:]-0.5)**2 + (y[:]-0.5)**2 - 0.5

        return displacement

    @mock.patch('coconut.tools.create_instance', side_effect=mock_create_instance)
    def test_displacement_on_nodes(self, create_instance):

        # test if nodes are moved to the correct position
        comb_solver = create_instance(self.comb_sol_par)
        model_part = comb_solver.get_interface_input().get_model_part(self.mp_name_in)
        x0, y0, z0 = model_part.x0, model_part.y0, model_part.z0
        displacement = self.get_displacement(x0, y0, z0)

        comb_solver.get_interface_input().set_variable_data(self.mp_name_in, 'displacement', displacement)

        # update position by iterating once in solver
        comb_solver.initialize()
        comb_solver.initialize_solution_step()
        comb_solver.solve_solution_step(comb_solver.get_interface_input())
        comb_solver.finalize_solution_step()
        comb_solver.finalize()

        master_disp = comb_solver.master_solver_wrapper.get_interface_input().get_interface_data()
        np.testing.assert_allclose(master_disp, displacement.ravel(), rtol=1e-12)

        for other_sol_wrapper in comb_solver.other_solver_wrapper_list:
            other_interface_input = other_sol_wrapper.get_interface_input()
            other_disp = other_interface_input.get_interface_data()
            np.testing.assert_allclose(other_disp, displacement.ravel(), rtol=1e-12)

    @mock.patch('coconut.tools.create_instance', side_effect=mock_create_instance)
    def test_output(self, create_instance):

        # test if nodes are moved to the correct position
        comb_solver = create_instance(self.comb_sol_par)
        model_part = comb_solver.get_interface_input().get_model_part(self.mp_name_in)
        x0, y0, z0 = model_part.x0, model_part.y0, model_part.z0
        displacement = self.get_displacement(x0, y0, z0)

        comb_solver.get_interface_input().set_variable_data(self.mp_name_in, 'displacement', displacement)

        # update position by iterating once in solver
        comb_solver.initialize()
        comb_solver.initialize_solution_step()
        output_interface = comb_solver.solve_solution_step(comb_solver.get_interface_input())
        comb_solver.finalize_solution_step()
        comb_solver.finalize()

        pressure = output_interface.get_variable_data(self.mp_name_out, 'pressure')
        traction = output_interface.get_variable_data(self.mp_name_out, 'traction')
        mp_out = comb_solver.get_interface_input().get_model_part(self.mp_name_out)
        disp_norm = np.linalg.norm(displacement.ravel())
        disp_min = np.min(displacement.ravel())
        pressure_value = disp_norm
        traction_value = [disp_min, -1 * disp_min, 2 * disp_min]
        pressure_ref = 2*np.full((mp_out.size,1), pressure_value)
        traction_ref = 2*np.full((mp_out.size, 3), traction_value)
        np.testing.assert_allclose(pressure, pressure_ref, rtol=1e-12)
        np.testing.assert_allclose(traction, traction_ref, rtol=1e-12)

if __name__ == '__main__':
    unittest.main()



