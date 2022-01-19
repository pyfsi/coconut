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
        object_module = __import__('coconut.tests.solver_wrappers.baffle.' + object_type, fromlist=[object_type])
    return object_module.create(settings)


class TestSolverWrapperBaffle(unittest.TestCase):

    def setUp(self):
        parameter_file_name = os.path.join(os.path.dirname(__file__), 'parameters.json')
        self.mul_factors= []

        with open(parameter_file_name, 'r') as parameter_file:
            parameters = json.load(parameter_file)

        self.baffle_sol_par = parameters

        self.mp_name_master_in = self.baffle_sol_par['settings']['master_model_part_inp']
        self.mp_name_slave_in = self.baffle_sol_par['settings']['slave_model_part_inp']
        self.mp_name_master_out = self.baffle_sol_par['settings']['master_model_part_out']
        self.mp_name_slave_out = self.baffle_sol_par['settings']['slave_model_part_out']


    def get_displacement(self, x, y, z):
        displacement = np.empty((x.size, 3))
        displacement[:, 0] = 0.0
        displacement[:, 1] = 0.0
        displacement[:, 2] = (x[:] - 0.5) ** 2 + (y[:] - 0.5) ** 2 - 0.5

        return displacement

    @mock.patch('coconut.tools.create_instance', side_effect=mock_create_instance)
    def test_displacement_on_nodes(self, create_instance):
        # test if nodes are moved to the correct position
        baffle_solver = create_instance(self.baffle_sol_par)
        baffle_solver.initialize()
        model_part = baffle_solver.get_interface_input().get_model_part(self.mp_name_master_in)
        x0, y0, z0 = model_part.x0, model_part.y0, model_part.z0
        displacement = self.get_displacement(x0, y0, z0)
        interface_input =  baffle_solver.get_interface_input()
        interface_input.set_variable_data(self.mp_name_master_in, 'displacement', displacement)

        # update position by iterating once in solver
        baffle_solver.initialize_solution_step()
        baffle_solver.solve_solution_step(baffle_solver.get_interface_input())
        baffle_solver.finalize_solution_step()
        baffle_solver.finalize()

        sol_wrapper_interface_input = baffle_solver.solver_wrapper.get_interface_input()
        master_disp = sol_wrapper_interface_input.get_variable_data(self.mp_name_master_in, 'displacement')
        slave_disp = sol_wrapper_interface_input.get_variable_data(self.mp_name_slave_in, 'displacement')
        np.testing.assert_allclose(master_disp, displacement, rtol=1e-12)
        np.testing.assert_allclose(slave_disp, displacement, rtol=1e-12)
        np.testing.assert_allclose(baffle_solver.get_interface_input().get_interface_data(), displacement.ravel(), rtol=1e-12)


    @mock.patch('coconut.tools.create_instance', side_effect=mock_create_instance)
    def test_output(self, create_instance):
        # test if nodes are moved to the correct position
        baffle_solver = create_instance(self.baffle_sol_par)
        baffle_solver.initialize()
        model_part = baffle_solver.get_interface_input().get_model_part(self.mp_name_master_in)
        x0, y0, z0 = model_part.x0, model_part.y0, model_part.z0
        displacement = self.get_displacement(x0, y0, z0)
        baffle_solver.get_interface_input().set_variable_data(self.mp_name_master_in, 'displacement', displacement)

        # update position by iterating once in solver
        baffle_solver.initialize_solution_step()
        output_interface = baffle_solver.solve_solution_step(baffle_solver.get_interface_input())
        baffle_solver.finalize_solution_step()
        baffle_solver.finalize()

        pressure = output_interface.get_variable_data(self.mp_name_master_out, 'pressure')
        traction = output_interface.get_variable_data(self.mp_name_master_out, 'traction')
        master_pres_ref, master_trac_ref = baffle_solver.solver_wrapper.calculate_output(displacement.ravel(),
                                                                                baffle_solver.get_interface_output(),
                                                                                self.mp_name_master_out)
        slave_pres_ref, slave_trac_ref = baffle_solver.solver_wrapper.calculate_output(displacement.ravel(),
                                                                                         baffle_solver.get_interface_output(),
                                                                                         self.mp_name_slave_out)

        pres_ref = master_pres_ref - slave_pres_ref
        trac_ref = master_trac_ref + slave_trac_ref

        np.testing.assert_allclose(pressure, pres_ref, rtol=1e-12)
        np.testing.assert_allclose(traction, trac_ref, rtol=1e-12)

if __name__ == '__main__':
    unittest.main()
