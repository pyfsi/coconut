from coconut import data_structure
import unittest
from coconut.data_structure.interface import Interface
from coconut.coupling_components.tools import create_instance
import json
import os
import numpy as np

class TestConvergenceCriterionAnd(unittest.TestCase):
    def test_convergence_criterion_and(self):
        m = 10
        dz = 2.0
        a0 = 1.0
        variable = 'area'
        mp_name = 'wall'
        interface_settings = [{'model_part': 'wall', 'variables': ['area']}]

        # Create interface

        model = data_structure.Model()
        ids = np.arange(0, m)
        x0 = np.zeros(m)
        y0 = np.zeros(m)
        z0 = np.arange(0, m * dz, dz)

        model.create_model_part(mp_name, x0, y0, z0, ids)

        a0_array = np.full((m, 1), a0)

        # create interface
        interface = Interface(interface_settings, model)
        interface.set_variable_data(mp_name, variable, a0_array)
        # model_part = model.CreateModelPart("wall")
        # model_part.AddNodalSolutionStepVariable(variable)
        # for i in range(m):
        #     model_part.CreateNewNode(i, 0.0, 0.0, i * dz)
        # step = 0
        # for node in model_part.Nodes:
        #     node.SetSolutionStepValue(variable, step, a0)
        # interface = Interface(model, interface_settings)

        parameter_file_name = os.path.join(os.path.dirname(__file__), 'test_and.json')
        with open(parameter_file_name, 'r') as parameter_file:
            settings = json.loads(parameter_file.read())
        convergence_criterion_and = create_instance(settings)
        convergence_criterion_and.initialize()
        for i in range(3):
            convergence_criterion_and.initialize_solution_step()
            is_satisfied = convergence_criterion_and.is_satisfied()
            self.assertFalse(is_satisfied)
            convergence_criterion_and.update(interface)
            is_satisfied = convergence_criterion_and.is_satisfied()
            self.assertFalse(is_satisfied)
            convergence_criterion_and.update(interface)
            is_satisfied = convergence_criterion_and.is_satisfied()
            self.assertFalse(is_satisfied)
            convergence_criterion_and.update(interface)
            is_satisfied = convergence_criterion_and.is_satisfied()
            self.assertFalse(is_satisfied)
            convergence_criterion_and.update(interface)
            is_satisfied = convergence_criterion_and.is_satisfied()
            self.assertTrue(is_satisfied)
            convergence_criterion_and.update(interface)
            is_satisfied = convergence_criterion_and.is_satisfied()
            self.assertTrue(is_satisfied)
            convergence_criterion_and.finalize_solution_step()


if __name__ == '__main__':
    unittest.main()
