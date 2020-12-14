from coconut import data_structure
import unittest
from coconut.data_structure.interface import Interface
from coconut.coupling_components.tools import create_instance

import os
import json
import numpy as np

class TestConvergenceCriterionRelativeNorm(unittest.TestCase):
    def test_convergence_criterion_relative_norm(self):
        m = 10
        dz = 2.0
        a0 = 10.0
        a1 = 1.0e-4
        a2 = 1.0e-6

        variable = 'area'
        model_part_name = 'wall'
        interface_settings = [{'model_part': 'wall', 'variables': ['area']}]

        # Create model and modelpart
        model = data_structure.Model()
        ids = np.arange(0, m)
        x0 = np.zeros(m)
        y0 = np.zeros(m)
        z0 = np.arange(0, m * dz, dz)
        mp_name = 'wall'

        model.create_model_part(mp_name, x0, y0, z0, ids)

        a0_array = np.full((m, 1), a0)
        a1_array = np.full((m, 1), a1)
        a2_array = np.full((m, 1), a2)
        # create interface
        interface = Interface(interface_settings, model)
        # interface_settings = json.loads('{"wall": "AREA"}')
        #
        # # Create interface
        # variable = vars(data_structure)["AREA"]
        # model = data_structure.Model()
        # model_part = model.CreateModelPart("wall")
        # model_part.AddNodalSolutionStepVariable(variable)
        # for i in range(m):
        #     model_part.CreateNewNode(i, 0.0, 0.0, i * dz)
        # step = 0
        # for node in model_part.Nodes:
        #     node.SetSolutionStepValue(variable, step, a0)
        # interface = Interface(model, interface_settings)

        parameter_file_name = os.path.join(os.path.dirname(__file__), 'test_relative_norm.json')
        with open(parameter_file_name, 'r') as parameter_file:
            settings = json.loads(parameter_file.read())

        convergence_criterion_relative_norm = create_instance(settings)
        convergence_criterion_relative_norm.Initialize()
        for i in range(3):
            convergence_criterion_relative_norm.InitializeSolutionStep()
            is_satisfied = convergence_criterion_relative_norm.IsSatisfied()
            self.assertFalse(is_satisfied)
            # for node in model_part.Nodes:
            #     node.SetSolutionStepValue(variable, step, a0)
            interface.set_variable_data(model_part_name, variable, a0_array)
            convergence_criterion_relative_norm.Update(interface)
            is_satisfied = convergence_criterion_relative_norm.IsSatisfied()
            self.assertFalse(is_satisfied)
            # for node in model_part.Nodes:
            #     node.SetSolutionStepValue(variable, step, a1)
            interface.set_variable_data(model_part_name, variable, a1_array)
            convergence_criterion_relative_norm.Update(interface)
            is_satisfied = convergence_criterion_relative_norm.IsSatisfied()
            self.assertFalse(is_satisfied)
            # for node in model_part.Nodes:
            #     node.SetSolutionStepValue(variable, step, a2)
            interface.set_variable_data(model_part_name, variable, a2_array)
            convergence_criterion_relative_norm.Update(interface)
            is_satisfied = convergence_criterion_relative_norm.IsSatisfied()
            self.assertTrue(is_satisfied)
            # for node in model_part.Nodes:
            #     node.SetSolutionStepValue(variable, step, a1)
            interface.set_variable_data(model_part_name, variable, a1_array)
            convergence_criterion_relative_norm.Update(interface)
            is_satisfied = convergence_criterion_relative_norm.IsSatisfied()
            self.assertFalse(is_satisfied)
            convergence_criterion_relative_norm.FinalizeSolutionStep()


if __name__ == '__main__':
    unittest.main()
