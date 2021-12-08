from coconut.tools import create_instance

import unittest
import os
import json
import numpy as np
from unittest import mock


def mock_create_instance(settings):
    object_type = settings['type']
    if 'dummy' not in object_type:
        return create_instance(settings)
    else:
        object_module = __import__('coconut.tests.solver_wrappers.explicit_combined.' + object_type, fromlist=[object_type])
    return object_module.create(settings)


class TestSolverWrapperExplicitCombined(unittest.TestCase):

    def setUp(self):
        parameter_file_name = os.path.join(os.path.dirname(__file__), 'parameters.json')
        self.mul_factors= []

        with open(parameter_file_name, 'r') as parameter_file:
            parameters = json.load(parameter_file)

        for sol_param in parameters['settings']['solver_wrappers']:
            if not sol_param['type'] == 'solver_wrappers.mapped':
                self.master_par_solver = sol_param
            else:
                self.mul_factors.append(sol_param['settings'].get('mul_factor', 1.0))

        self.exp_comb_sol_par = parameters

        self.mp_name_in = self.master_par_solver['settings']['interface_input'][0]['model_part']
        self.mp_name_out = self.master_par_solver['settings']['interface_output'][0]['model_part']

    def get_displacement(self, x, y, z):
        displacement = np.empty((x.size, 3))
        displacement[:, 0] = 0.0
        displacement[:, 1] = 0.0
        displacement[:, 2] = (x[:] - 0.5) ** 2 + (y[:] - 0.5) ** 2 - 0.5

        return displacement

    @mock.patch('coconut.tools.create_instance', side_effect=mock_create_instance)
    def test_displacement_on_nodes(self, create_instance):
        # test if nodes are moved to the correct position
        exp_comb_solver = create_instance(self.exp_comb_sol_par)
        model_part = exp_comb_solver.get_interface_input().get_model_part(self.mp_name_in)
        x0, y0, z0 = model_part.x0, model_part.y0, model_part.z0
        displacement_1 = self.get_displacement(x0, y0, z0)

        exp_comb_solver.get_interface_input().set_variable_data(self.mp_name_in, 'displacement', displacement_1)

        # update position by iterating once in solver
        exp_comb_solver.initialize()
        exp_comb_solver.initialize_solution_step()
        # Solving for 1st iteration
        exp_comb_solver.solve_solution_step(exp_comb_solver.get_interface_input())
        interface_input = exp_comb_solver.master_solver_wrapper.get_interface_input()
        out_disp = interface_input.get_interface_data()
        np.testing.assert_allclose(out_disp, displacement_1.ravel(), rtol=1e-12)
        for index, mapped_sol_wrapper in enumerate(exp_comb_solver.solver_wrapper_list):
            if not index == exp_comb_solver.master_sol_index:
                interface_input = mapped_sol_wrapper.solver_wrapper.get_interface_input()
                out_disp = interface_input.get_interface_data()
                np.testing.assert_allclose(out_disp, displacement_1.ravel(), rtol=1e-12)

        # Solving for 2st iteration
        displacement_2 = 2*displacement_1
        exp_comb_solver.get_interface_input().set_variable_data(self.mp_name_in, 'displacement', displacement_2)
        exp_comb_solver.solve_solution_step(exp_comb_solver.get_interface_input())
        interface_input = exp_comb_solver.master_solver_wrapper.get_interface_input()
        out_disp = interface_input.get_interface_data()
        np.testing.assert_allclose(out_disp, displacement_2.ravel(), rtol=1e-12)
        for index, mapped_sol_wrapper in enumerate(exp_comb_solver.solver_wrapper_list):
            if not index == exp_comb_solver.master_sol_index:
                interface_input = mapped_sol_wrapper.solver_wrapper.get_interface_input()
                out_disp = interface_input.get_interface_data()
                # mapped interfaces are not updated
                np.testing.assert_allclose(out_disp, displacement_1.ravel(), rtol=1e-12)

        exp_comb_solver.finalize_solution_step()
        exp_comb_solver.finalize()

    @mock.patch('coconut.tools.create_instance', side_effect=mock_create_instance)
    def test_output(self, create_instance):

        exp_comb_solver = create_instance(self.exp_comb_sol_par)
        model_part = exp_comb_solver.get_interface_input().get_model_part(self.mp_name_in)
        x0, y0, z0 = model_part.x0, model_part.y0, model_part.z0
        displacement_1 = self.get_displacement(x0, y0, z0)
        exp_comb_solver.get_interface_input().set_variable_data(self.mp_name_in, 'displacement', displacement_1)
        # update position by iterating once in solver
        exp_comb_solver.initialize()
        exp_comb_solver.initialize_solution_step()
        # Solving for 1st iteration
        output_interface = exp_comb_solver.solve_solution_step(exp_comb_solver.get_interface_input())
        pressure_1 = output_interface.get_variable_data(self.mp_name_out, 'pressure')
        traction_1 = output_interface.get_variable_data(self.mp_name_out, 'traction')
        pres_ref_1, trac_ref_1 = exp_comb_solver.master_solver_wrapper.calculate_output(displacement_1.ravel(),
                                                                                    exp_comb_solver.get_interface_output(),
                                                                                self.mp_name_out)
        for mul_factor, mapped_sol_wrapper in zip(self.mul_factors,exp_comb_solver.mapped_solver_wrapper_list):
            pres_other_1, trac_other_1 = mapped_sol_wrapper.solver_wrapper.calculate_output(displacement_1.ravel(),
                                                                                        exp_comb_solver.get_interface_output(),
                                                                                        self.mp_name_out)
            pres_ref_1 += mul_factor*pres_other_1
            trac_ref_1 += mul_factor*trac_other_1

        np.testing.assert_allclose(pressure_1, pres_ref_1, rtol=1e-12)
        np.testing.assert_allclose(traction_1, trac_ref_1, rtol=1e-12)

        # Solving for 2nd iteration
        displacement_2 = 2*displacement_1
        exp_comb_solver.get_interface_input().set_variable_data(self.mp_name_in, 'displacement', displacement_2)
        output_interface = exp_comb_solver.solve_solution_step(exp_comb_solver.get_interface_input())
        pressure_2 = output_interface.get_variable_data(self.mp_name_out, 'pressure')
        traction_2 = output_interface.get_variable_data(self.mp_name_out, 'traction')
        pres_ref_2, trac_ref_2 = exp_comb_solver.master_solver_wrapper.calculate_output(displacement_2.ravel(),
                                                                                    exp_comb_solver.get_interface_output(),
                                                                                    self.mp_name_out)
        for mul_factor, mapped_sol_wrapper in zip(self.mul_factors,exp_comb_solver.mapped_solver_wrapper_list):
            #only first iteration values of iterface_other are added
            pres_ref_2 += mul_factor*pres_other_1
            trac_ref_2 += mul_factor*trac_other_1

        np.testing.assert_allclose(pressure_2, pres_ref_2, rtol=1e-12)
        np.testing.assert_allclose(traction_2, trac_ref_2, rtol=1e-12)

        exp_comb_solver.finalize_solution_step()
        exp_comb_solver.finalize()

if __name__ == '__main__':
    unittest.main()
