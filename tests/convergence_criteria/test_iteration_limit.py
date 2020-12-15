from coconut import data_structure
from coconut.data_structure.interface import Interface
from coconut.coupling_components.tools import create_instance

import unittest
import os
import json
import numpy as np


class TestConvergenceCriterionIterationLimit(unittest.TestCase):

    def test_convergence_criterion_iteration_limit(self):
        m = 10
        dz = 2
        a0 = 1
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

        # create interface
        interface = Interface(interface_settings, model)
        interface.set_variable_data(model_part_name, variable, a0_array)

        # read settings
        parameter_file_name = os.path.join(os.path.dirname(__file__), 'test_iteration_limit.json')
        with open(parameter_file_name, 'r') as parameter_file:
            settings = json.load(parameter_file)

        convergence_criterion_iteration_limit = create_instance(settings)
        convergence_criterion_iteration_limit.initialize()
        for i in range(3):
            convergence_criterion_iteration_limit.initialize_solution_step()
            is_satisfied = convergence_criterion_iteration_limit.is_satisfied()
            self.assertFalse(is_satisfied)
            convergence_criterion_iteration_limit.update(interface)
            is_satisfied = convergence_criterion_iteration_limit.is_satisfied()
            self.assertFalse(is_satisfied)
            convergence_criterion_iteration_limit.update(interface)
            is_satisfied = convergence_criterion_iteration_limit.is_satisfied()
            self.assertTrue(is_satisfied)
            convergence_criterion_iteration_limit.update(interface)
            is_satisfied = convergence_criterion_iteration_limit.is_satisfied()
            self.assertTrue(is_satisfied)
            convergence_criterion_iteration_limit.finalize_solution_step()


if __name__ == '__main__':
    unittest.main()
