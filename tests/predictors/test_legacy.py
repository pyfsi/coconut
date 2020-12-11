from coconut import data_structure
import unittest
from coconut.coupling_components.tools import CreateInstance
from coconut.coupling_components.interface import Interface

import numpy as np


class TestPredictorLegacy(unittest.TestCase):
    def test_predictor_legacy(self):
        m = 10
        dz = 3.0
        a0 = 5.0
        p1 = 5.0
        a1 = 1.0
        p2 = -3.0
        a2 = 5.0
        p3 = 13.0
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

        # Create predictor
        parameter_file_name = "predictors/test_legacy.json"
        with open(parameter_file_name, 'r') as parameter_file:
            settings = data_structure.Parameters(parameter_file.read())

        predictor_legacy = CreateInstance(settings)
        predictor_legacy.Initialize(interface)

        # Test predictor: first prediction needs to be equal to initialized value
        predictor_legacy.InitializeSolutionStep()
        prediction = predictor_legacy.Predict(interface)
        self.assertIsInstance(prediction, Interface)
        prediction_as_array = prediction.GetNumpyArray()
        for i in range(m):
            self.assertAlmostEqual(p1, prediction_as_array[i])
        interface_as_array = a1 * np.ones_like(prediction_as_array)
        interface.SetNumpyArray(interface_as_array)
        predictor_legacy.Update(interface)
        predictor_legacy.FinalizeSolutionStep()

        # Test predictor: second prediction needs to be linear
        predictor_legacy.InitializeSolutionStep()
        prediction = predictor_legacy.Predict(interface)
        self.assertIsInstance(prediction, Interface)
        prediction_as_array = prediction.GetNumpyArray()
        for i in range(m):
            self.assertAlmostEqual(p2, prediction_as_array[i])
        interface_as_array = a2 * np.ones_like(prediction_as_array)
        interface.SetNumpyArray(interface_as_array)
        predictor_legacy.Update(interface)
        predictor_legacy.FinalizeSolutionStep()

        # Test predictor: third prediction needs to be legacy
        predictor_legacy.InitializeSolutionStep()
        prediction = predictor_legacy.Predict(interface)
        self.assertIsInstance(prediction, Interface)
        prediction_as_array = prediction.GetNumpyArray()
        for i in range(m):
            self.assertAlmostEqual(p3, prediction_as_array[i])


if __name__ == '__main__':
    unittest.main()
