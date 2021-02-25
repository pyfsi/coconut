from coconut import data_structure
from coconut.data_structure.interface import Interface
from coconut.tools import create_instance

import unittest
import numpy as np


class TestConvergenceCriterionRelativeNorm(unittest.TestCase):

    def test_convergence_criterion_relative_norm(self):
        m = 10
        dz = 2
        a0 = 10
        a1 = 1e-4
        a2 = 1e-6
        variable = 'area'
        model_part_name = 'wall'
        interface_settings = [{'model_part': model_part_name, 'variables': [variable]}]

        # create model and model_part
        model = data_structure.Model()
        ids = np.arange(0, m)
        x0 = np.zeros(m)
        y0 = np.zeros(m)
        z0 = np.arange(0, m * dz, dz)
        model.create_model_part(model_part_name, x0, y0, z0, ids)

        a0_array = np.full((m, 1), a0)
        a1_array = np.full((m, 1), a1)
        a2_array = np.full((m, 1), a2)

        # create interface
        interface = Interface(interface_settings, model)

        # create convergence criterion
        parameters = {'type': 'convergence_criteria.relative_norm',
                      'settings': {'tolerance': 1e-6, 'order': 2}}
        convergence_criterion_relative_norm = create_instance(parameters)
        convergence_criterion_relative_norm.initialize()

        # test convergence criterion
        for i in range(3):
            convergence_criterion_relative_norm.initialize_solution_step()
            is_satisfied = convergence_criterion_relative_norm.is_satisfied()
            self.assertFalse(is_satisfied)
            interface.set_variable_data(model_part_name, variable, a0_array)
            convergence_criterion_relative_norm.update(interface)
            is_satisfied = convergence_criterion_relative_norm.is_satisfied()
            self.assertFalse(is_satisfied)
            interface.set_variable_data(model_part_name, variable, a1_array)
            convergence_criterion_relative_norm.update(interface)
            is_satisfied = convergence_criterion_relative_norm.is_satisfied()
            self.assertFalse(is_satisfied)
            interface.set_variable_data(model_part_name, variable, a2_array)
            convergence_criterion_relative_norm.update(interface)
            is_satisfied = convergence_criterion_relative_norm.is_satisfied()
            self.assertTrue(is_satisfied)
            interface.set_variable_data(model_part_name, variable, a1_array)
            convergence_criterion_relative_norm.update(interface)
            is_satisfied = convergence_criterion_relative_norm.is_satisfied()
            self.assertFalse(is_satisfied)
            convergence_criterion_relative_norm.finalize_solution_step()


if __name__ == '__main__':
    unittest.main()
