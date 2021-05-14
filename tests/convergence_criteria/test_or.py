from coconut import data_structure
from coconut.data_structure.interface import Interface
from coconut.tools import create_instance

import unittest
import json
import os
import numpy as np


class TestConvergenceCriterionOr(unittest.TestCase):

    def test_convergence_criterion_or(self):
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
        parameter_file_name = os.path.join(os.path.dirname(__file__), 'test_or.json')
        with open(parameter_file_name) as parameter_file:
            parameters = json.load(parameter_file)

        convergence_criterion_or = create_instance(parameters)
        convergence_criterion_or.initialize()
        for i in range(3):
            convergence_criterion_or.initialize_solution_step()
            is_satisfied = convergence_criterion_or.is_satisfied()
            self.assertFalse(is_satisfied)
            convergence_criterion_or.update(interface)
            is_satisfied = convergence_criterion_or.is_satisfied()
            self.assertFalse(is_satisfied)
            convergence_criterion_or.update(interface)
            is_satisfied = convergence_criterion_or.is_satisfied()
            self.assertTrue(is_satisfied)
            convergence_criterion_or.update(interface)
            is_satisfied = convergence_criterion_or.is_satisfied()
            self.assertTrue(is_satisfied)
            convergence_criterion_or.update(interface)
            is_satisfied = convergence_criterion_or.is_satisfied()
            self.assertTrue(is_satisfied)
            convergence_criterion_or.update(interface)
            is_satisfied = convergence_criterion_or.is_satisfied()
            self.assertTrue(is_satisfied)
            convergence_criterion_or.finalize_solution_step()


if __name__ == '__main__':
    unittest.main()
