from coconut import data_structure
import unittest
from coconut.coupling_components.tools import create_instance
from coconut.data_structure.interface import Interface

import os
import numpy as np
import json

class TestPredictorQuadratic(unittest.TestCase):
    def test_predictor_quadratic(self):
        m = 10
        dz = 3.0
        a0 = 5.0
        p1 = 5.0
        a1 = 1.0
        p2 = -3.0
        a2 = 5.0
        p3 = 17.0
        # interface_settings = data_structure.Parameter('{"wall": "AREA"}')
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

        # Create predictor
        parameter_file_name = os.path.join(os.path.dirname(__file__),"test_quadratic.json")
        with open(parameter_file_name, 'r') as parameter_file:
            settings = json.loads(parameter_file.read())

        predictor_quadratic = create_instance(settings)
        predictor_quadratic.initialize(interface)

        # Test predictor: first prediction needs to be equal to initialized value
        predictor_quadratic.initialize_solution_step()
        prediction = predictor_quadratic.predict(interface)
        self.assertIsInstance(prediction, Interface)
        prediction_as_array = prediction.get_interface_data()
        for i in range(m):
            self.assertAlmostEqual(p1, prediction_as_array[i])
        interface_as_array = a1 * np.ones_like(prediction_as_array)
        interface.set_interface_data(interface_as_array)
        predictor_quadratic.update(interface)
        predictor_quadratic.finalize_solution_step()

        # Test predictor: second prediction needs to be linear
        predictor_quadratic.initialize_solution_step()
        prediction = predictor_quadratic.predict(interface)
        self.assertIsInstance(prediction, Interface)
        prediction_as_array = prediction.get_interface_data()
        for i in range(m):
            self.assertAlmostEqual(p2, prediction_as_array[i])
        interface_as_array = a2 * np.ones_like(prediction_as_array)
        interface.set_interface_data(interface_as_array)
        predictor_quadratic.update(interface)
        predictor_quadratic.finalize_solution_step()

        # Test predictor: third prediction needs to be quadratic
        predictor_quadratic.initialize_solution_step()
        prediction = predictor_quadratic.predict(interface)
        self.assertIsInstance(prediction, Interface)
        prediction_as_array = prediction.get_interface_data()
        for i in range(m):
            self.assertAlmostEqual(p3, prediction_as_array[i])


if __name__ == '__main__':
    unittest.main()
