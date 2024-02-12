from coconut.tools import create_instance

import unittest
import os
import json
import numpy as np
from unittest import mock

from coconut import tools


def mock_create_instance(settings):
    object_type = settings['type']
    if 'dummy' not in object_type:
        return create_instance(settings)
    else:
        object_module = __import__('coconut.tests.solver_wrappers.combined.' + object_type, fromlist=[object_type])
    return object_module.create(settings)


class TestSolverWrapperCombined(unittest.TestCase):

    def setUp(self):
        parameter_file_name = os.path.join(os.path.dirname(__file__), 'parameters.json')

        with open(parameter_file_name, 'r') as parameter_file:
            parameters = json.load(parameter_file)

        for sol_param in parameters['settings']['solver_wrappers']:
            if not sol_param["type"] == "solver_wrappers.mapped":
                self.master_par_solver = sol_param

        self.comb_sol_par = parameters

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
        comb_solver = create_instance(self.comb_sol_par)
        model_part = comb_solver.get_interface_input().get_model_part(self.mp_name_in)
        x0, y0, z0 = model_part.x0, model_part.y0, model_part.z0
        displacement = self.get_displacement(x0, y0, z0)

        interface_input = comb_solver.get_interface_input()
        interface_input.set_variable_data(self.mp_name_in, 'displacement', displacement)

        # update position by iterating once in solver
        comb_solver.initialize()
        comb_solver.initialize_solution_step()
        comb_solver.solve_solution_step(interface_input)
        comb_solver.finalize_solution_step()
        comb_solver.finalize()

        interface_input = comb_solver.master_solver_wrapper.get_interface_input()
        out_disp = interface_input.get_interface_data()
        np.testing.assert_allclose(out_disp, displacement.ravel(), rtol=1e-12)

        for index, mapped_sol_wrapper in enumerate(comb_solver.mapped_solver_wrappers):
            interface_input = mapped_sol_wrapper.solver_wrapper.get_interface_input()
            out_disp = interface_input.get_interface_data()
            np.testing.assert_allclose(out_disp, displacement.ravel(), rtol=1e-12)

    @mock.patch('coconut.tools.create_instance', side_effect=mock_create_instance)
    def test_output(self, create_instance):
        # test if nodes are moved to the correct position
        comb_solver = create_instance(self.comb_sol_par)
        model_part = comb_solver.get_interface_input().get_model_part(self.mp_name_in)
        x0, y0, z0 = model_part.x0, model_part.y0, model_part.z0
        displacement = self.get_displacement(x0, y0, z0)

        interface_input = comb_solver.get_interface_input()
        interface_input.set_variable_data(self.mp_name_in, 'displacement', displacement)

        # update position by iterating once in solver
        comb_solver.initialize()
        comb_solver.initialize_solution_step()
        output_interface = comb_solver.solve_solution_step(interface_input)
        comb_solver.finalize_solution_step()
        comb_solver.finalize()

        pressure = output_interface.get_variable_data(self.mp_name_out, 'pressure')
        traction = output_interface.get_variable_data(self.mp_name_out, 'traction')
        pres_ref, trac_ref = \
            comb_solver.master_solver_wrapper.calculate_output(displacement.ravel(), comb_solver.get_interface_output(),
                                                               self.mp_name_out)
        for index, mapped_sol_wrapper in enumerate(comb_solver.mapped_solver_wrappers):
            pres_other, trac_other = \
                mapped_sol_wrapper.solver_wrapper.calculate_output(displacement.ravel(),
                                                                   comb_solver.get_interface_output(), self.mp_name_out)
            pres_ref += pres_other
            trac_ref += trac_other

        np.testing.assert_allclose(pressure, pres_ref, rtol=1e-12)
        np.testing.assert_allclose(traction, trac_ref, rtol=1e-12)

        tools.print_info('\n')
        comb_solver.print_components_info('├─')
        tools.print_info('\n')
        tools.print_info(comb_solver.get_time_allocation())


if __name__ == '__main__':
    unittest.main()
