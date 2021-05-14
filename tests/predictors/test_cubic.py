from coconut import data_structure
from coconut.data_structure.interface import Interface
from coconut.tools import create_instance

import unittest
import os
import json
import numpy as np


class TestPredictorCubic(unittest.TestCase):

    def test_predictor_cubic(self):
        m = 10
        dz = 3
        a0 = 1
        p1 = 1
        a1 = 3
        p2 = 5
        a2 = 37
        p3 = 103
        a3 = 151
        p4 = 393
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

        # create predictor
        parameters = {'type': 'predictors.cubic'}
        predictor_cubic = create_instance(parameters)
        predictor_cubic.initialize(interface)

        # first prediction needs to be equal to initialized value
        predictor_cubic.initialize_solution_step()
        prediction = predictor_cubic.predict(interface)
        self.assertIsInstance(prediction, Interface)
        prediction_as_array = prediction.get_interface_data()
        for i in range(m):
            self.assertAlmostEqual(p1, prediction_as_array[i])
        interface_as_array = a1 * np.ones_like(prediction_as_array)
        interface.set_interface_data(interface_as_array)
        predictor_cubic.update(interface)
        predictor_cubic.finalize_solution_step()

        # second prediction needs to be linear
        predictor_cubic.initialize_solution_step()
        prediction = predictor_cubic.predict(interface)
        self.assertIsInstance(prediction, Interface)
        prediction_as_array = prediction.get_interface_data()
        for i in range(m):
            self.assertAlmostEqual(p2, prediction_as_array[i])
        interface_as_array = a2 * np.ones_like(prediction_as_array)
        interface.set_interface_data(interface_as_array)
        predictor_cubic.update(interface)
        predictor_cubic.finalize_solution_step()

        # third prediction needs to be quadratic
        predictor_cubic.initialize_solution_step()
        prediction = predictor_cubic.predict(interface)
        self.assertIsInstance(prediction, Interface)
        prediction_as_array = prediction.get_interface_data()
        for i in range(m):
            self.assertAlmostEqual(p3, prediction_as_array[i])
        interface_as_array = a3 * np.ones_like(prediction_as_array)
        interface.set_interface_data(interface_as_array)
        predictor_cubic.update(interface)
        predictor_cubic.finalize_solution_step()
        
        # fourth prediction needs to be cubic
        predictor_cubic.initialize_solution_step()
        prediction = predictor_cubic.predict(interface)
        self.assertIsInstance(prediction, Interface)
        prediction_as_array = prediction.get_interface_data()
        for i in range(m):
            self.assertAlmostEqual(p4, prediction_as_array[i])


if __name__ == '__main__':
    unittest.main()
