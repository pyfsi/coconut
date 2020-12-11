from coconut import data_structure
import unittest
from coconut.coupling_components.interface import Interface
from coconut.coupling_components.tools import CreateInstance

import os


class TestConvergenceCriterionAbsoluteNorm(unittest.TestCase):
    def test_convergence_criterion_absolute_norm(self):
        m = 10
        dz = 2.0
        a0 = 10.0
        a1 = 1.0e-4
        a2 = 1.0e-7
        interface_settings = data_structure.Parameters('{"wall": "AREA"}')

        # Create interface
        variable = vars(data_structure)["AREA"]
        model = data_structure.Model()
        model_part = model.CreateModelPart("wall")
        model_part.AddNodalSolutionStepVariable(variable)
        for i in range(m):
            model_part.CreateNewNode(i, 0.0, 0.0, i * dz)
        step = 0
        for node in model_part.Nodes:
            node.SetSolutionStepValue(variable, step, a0)
        interface = Interface(model, interface_settings)

        parameter_file_name = os.path.join(os.path.dirname(__file__), 'test_absolute_norm.json')
        with open(parameter_file_name, 'r') as parameter_file:
            settings = data_structure.Parameters(parameter_file.read())

        convergence_criterion_absolute_norm = CreateInstance(settings)
        convergence_criterion_absolute_norm.Initialize()
        for i in range(3):
            convergence_criterion_absolute_norm.InitializeSolutionStep()
            is_satisfied = convergence_criterion_absolute_norm.IsSatisfied()
            self.assertFalse(is_satisfied)
            for node in model_part.Nodes:
                node.SetSolutionStepValue(variable, step, a0)
            convergence_criterion_absolute_norm.Update(interface)
            is_satisfied = convergence_criterion_absolute_norm.IsSatisfied()
            self.assertFalse(is_satisfied)
            for node in model_part.Nodes:
                node.SetSolutionStepValue(variable, step, a1)
            convergence_criterion_absolute_norm.Update(interface)
            is_satisfied = convergence_criterion_absolute_norm.IsSatisfied()
            self.assertFalse(is_satisfied)
            for node in model_part.Nodes:
                node.SetSolutionStepValue(variable, step, a2)
            convergence_criterion_absolute_norm.Update(interface)
            is_satisfied = convergence_criterion_absolute_norm.IsSatisfied()
            self.assertTrue(is_satisfied)
            for node in model_part.Nodes:
                node.SetSolutionStepValue(variable, step, a1)
            convergence_criterion_absolute_norm.Update(interface)
            is_satisfied = convergence_criterion_absolute_norm.IsSatisfied()
            self.assertFalse(is_satisfied)
            convergence_criterion_absolute_norm.FinalizeSolutionStep()


if __name__ == '__main__':
    unittest.main()
