from coconut import data_structure
import unittest
from coconut.data_structure.interface import Interface
from coconut.coupling_components.tools import create_instance

import os
import json
import numpy as np

class TestConvergenceCriterionIterationLimit(unittest.TestCase):
    def test_convergence_criterion_iteration_limit(self):
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
        # interface_settings = json.loads('{"wall": "AREA"}')
        #
        # # Create interface
        # variable = vars(data_structure)["AREA"]
        # model = data_structure.Model()
        # model_part = model.CreateModelPart("wall")
        # model_part.AddNodalSolutionStepVariable(variable)
        # for i in range(m):
        #     model_part.CreateNewNode(i, 0.0, 0.0, i*dz)
        # step = 0
        # for node in model_part.Nodes:
        #     node.SetSolutionStepValue(variable, step, a0)
        # interface = Interface(model, interface_settings)

        parameter_file_name = os.path.join(os.path.dirname(__file__), 'test_iteration_limit.json')
        with open(parameter_file_name, 'r') as parameter_file:
            settings = json.loads(parameter_file.read())

        convergence_criterion_iteration_limit = create_instance(settings)
        convergence_criterion_iteration_limit.Initialize()
        for i in range(3):
            convergence_criterion_iteration_limit.InitializeSolutionStep()
            is_satisfied = convergence_criterion_iteration_limit.IsSatisfied()
            self.assertFalse(is_satisfied)
            convergence_criterion_iteration_limit.Update(interface)
            is_satisfied = convergence_criterion_iteration_limit.IsSatisfied()
            self.assertFalse(is_satisfied)
            convergence_criterion_iteration_limit.Update(interface)
            is_satisfied = convergence_criterion_iteration_limit.IsSatisfied()
            self.assertTrue(is_satisfied)
            convergence_criterion_iteration_limit.Update(interface)
            is_satisfied = convergence_criterion_iteration_limit.IsSatisfied()
            self.assertTrue(is_satisfied)
            convergence_criterion_iteration_limit.FinalizeSolutionStep()


if __name__ == '__main__':
    unittest.main()
