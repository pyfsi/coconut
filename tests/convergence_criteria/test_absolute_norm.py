from coconut import data_structure
import unittest
from coconut.data_structure.interface import Interface
from coconut.coupling_components.tools import create_instance

import os
import json
import numpy as np

class TestConvergenceCriterionAbsoluteNorm(unittest.TestCase):
    def test_convergence_criterion_absolute_norm(self):
        m = 10
        dz = 2.0
        a0 = 10.0
        a1 = 1.0e-4
        a2 = 1.0e-7
        variable = 'area'
        model_part_name = 'wall'
        interface_settings = [{'model_part': 'wall', 'variables': ['area']}]


        # Create model and modelpart
        model = data_structure.Model()
        ids = np.arange(0,m)
        x0 = np.zeros(m)
        y0 = np.zeros(m)
        z0 = np.arange(0,m * dz,dz)
        mp_name = 'wall'

        model.create_model_part(mp_name, x0, y0, z0, ids)

        a0_array = np.full((m,1), a0)
        a1_array = np.full((m, 1), a1)
        a2_array = np.full((m, 1), a2)
        # create interface
        interface = Interface(interface_settings, model)
        # model_part.AddNodalSolutionStepVariable(variable)
        # for i in range(m):
        #     model_part.CreateNewNode(i, 0.0, 0.0, i * dz)
        # step = 0
        # for node in model_part.Nodes:
        #     node.SetSolutionStepValue(variable, step, a0)
        # interface = Interface(model, interface_settings)

        parameter_file_name = os.path.join(os.path.dirname(__file__), 'test_absolute_norm.json')
        with open(parameter_file_name, 'r') as parameter_file:
            settings = json.loads(parameter_file.read())

        convergence_criterion_absolute_norm = create_instance(settings)
        convergence_criterion_absolute_norm.initialize()
        for i in range(3):
            convergence_criterion_absolute_norm.initialize_solution_step()
            is_satisfied = convergence_criterion_absolute_norm.is_satisfied()
            self.assertFalse(is_satisfied)
            # for node in model_part.Nodes:
            #     node.SetSolutionStepValue(variable, step, a0)
            interface.set_variable_data(model_part_name, variable, a0_array)
            convergence_criterion_absolute_norm.update(interface)
            is_satisfied = convergence_criterion_absolute_norm.is_satisfied()
            self.assertFalse(is_satisfied)
            # for node in model_part.Nodes:
            #     node.SetSolutionStepValue(variable, step, a1)
            interface.set_variable_data(model_part_name, variable, a1_array)
            convergence_criterion_absolute_norm.update(interface)
            is_satisfied = convergence_criterion_absolute_norm.is_satisfied()
            self.assertFalse(is_satisfied)
            # for node in model_part.Nodes:
            #     node.SetSolutionStepValue(variable, step, a2)
            interface.set_variable_data(model_part_name, variable, a2_array)
            convergence_criterion_absolute_norm.update(interface)
            is_satisfied = convergence_criterion_absolute_norm.is_satisfied()
            self.assertTrue(is_satisfied)
            # for node in model_part.Nodes:
            #     node.SetSolutionStepValue(variable, step, a1)
            interface.set_variable_data(model_part_name, variable, a1_array)
            convergence_criterion_absolute_norm.update(interface)
            is_satisfied = convergence_criterion_absolute_norm.is_satisfied()
            self.assertFalse(is_satisfied)
            convergence_criterion_absolute_norm.finalize_solution_step()


if __name__ == '__main__':
    unittest.main()
