from coconut import data_structure
import unittest
from coconut.coupling_components.tools import create_instance
from coconut.data_structure.interface import Interface

import os
import numpy as np
import json
class TestPredictor(unittest.TestCase):
    def test_predictor(self):
        m = 10
        dz = 3.0
        a0 = 1.0
        a1 = 2.0
        a2 = 3.0
        a3 = 4.0
        a4 = 5.0

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
        #     model_part.CreateNewNode(i, 0.0, 0.0, i * dz)
        # step = 0
        # for node in model_part.Nodes:
        #     node.SetSolutionStepValue(variable, step, a0)
        # interface = Interface(model, interface_settings)

        # Create predictor
        parameter_file_name = os.path.join(os.path.dirname(__file__ ),"test_cubic.json")
        with open(parameter_file_name, 'r') as parameter_file:
            settings = json.loads(parameter_file.read())
        predictor_cubic = create_instance(settings)
        predictor_cubic.Initialize(interface)
        interface_as_array = interface.get_interface_data()

        # Test predictor: a linear relation should be predicted in the same way
        # by linear, quadratic and cubic predictors
        predictor_cubic.InitializeSolutionStep()
        interface.set_interface_data(a1 * interface_as_array)
        predictor_cubic.Update(interface)
        predictor_cubic.FinalizeSolutionStep()
        predictor_cubic.InitializeSolutionStep()
        interface.set_interface_data(a2 * interface_as_array)
        predictor_cubic.Update(interface)
        predictor_cubic.FinalizeSolutionStep()
        predictor_cubic.InitializeSolutionStep()
        interface.set_interface_data(a3 * interface_as_array)
        predictor_cubic.Update(interface)
        predictor_cubic.FinalizeSolutionStep()

        predictor_cubic.InitializeSolutionStep()
        prediction_linear = predictor_cubic.linear(interface).get_interface_data()
        prediction_quadratic = predictor_cubic.quadratic(interface).get_interface_data()
        prediction_cubic = predictor_cubic.cubic(interface).get_interface_data()
        for i in range(m):
            self.assertAlmostEqual(a4, prediction_linear[i])
            self.assertAlmostEqual(a4, prediction_quadratic[i])
            self.assertAlmostEqual(a4, prediction_cubic[i])

        # Test predictor: error if no update
        with self.assertRaises(Exception):
            predictor_cubic.InitializeSolutionStep()
            predictor_cubic.FinalizeSolutionStep()

        # Test predictor: error if updated twice
        with self.assertRaises(Exception):
            predictor_cubic.InitializeSolutionStep()
            prediction = predictor_cubic.predict(interface)
            prediction = predictor_cubic.predict(interface)
            predictor_cubic.FinalizeSolutionStep()

        # Test predictor: error if prediction after update
        with self.assertRaises(Exception):
            predictor_cubic.InitializeSolutionStep()
            prediction = predictor_cubic.Update(interface)
            prediction = predictor_cubic.predict(interface)
            predictor_cubic.FinalizeSolutionStep()


if __name__ == '__main__':
    unittest.main()
