from coconut import data_structure
from coconut.data_structure.interface import Interface
from coconut.tools import create_instance

import unittest
import numpy as np


class TestPredictor(unittest.TestCase):

    def test_predictor(self):
        m = 10
        dz = 3
        a0 = 1
        a1 = 2
        a2 = 3
        a3 = 4
        a4 = 5
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
        interface_as_array = interface.get_interface_data()

        # create predictor
        parameters = {'type': 'predictors.cubic'}
        predictor_cubic = create_instance(parameters)
        predictor_cubic.initialize(interface)

        # a linear relation should be predicted in the same way by linear, quadratic and cubic predictors
        predictor_cubic.initialize_solution_step()
        interface.set_interface_data(a1 * interface_as_array)
        predictor_cubic.update(interface)
        predictor_cubic.finalize_solution_step()
        predictor_cubic.initialize_solution_step()
        interface.set_interface_data(a2 * interface_as_array)
        predictor_cubic.update(interface)
        predictor_cubic.finalize_solution_step()
        predictor_cubic.initialize_solution_step()
        interface.set_interface_data(a3 * interface_as_array)
        predictor_cubic.update(interface)
        predictor_cubic.finalize_solution_step()

        predictor_cubic.initialize_solution_step()
        prediction_linear = predictor_cubic.linear(interface).get_interface_data()
        prediction_quadratic = predictor_cubic.quadratic(interface).get_interface_data()
        prediction_cubic = predictor_cubic.cubic(interface).get_interface_data()
        for i in range(m):
            self.assertAlmostEqual(a4, prediction_linear[i])
            self.assertAlmostEqual(a4, prediction_quadratic[i])
            self.assertAlmostEqual(a4, prediction_cubic[i])

        # rror if no update
        with self.assertRaises(Exception):
            predictor_cubic.initialize_solution_step()
            predictor_cubic.finalize_solution_step()

        # error if updated twice
        with self.assertRaises(Exception):
            predictor_cubic.initialize_solution_step()
            _ = predictor_cubic.predict(interface)
            _ = predictor_cubic.predict(interface)
            predictor_cubic.finalize_solution_step()

        # error if prediction after update
        with self.assertRaises(Exception):
            predictor_cubic.initialize_solution_step()
            _ = predictor_cubic.update(interface)
            _ = predictor_cubic.predict(interface)
            predictor_cubic.finalize_solution_step()


if __name__ == '__main__':
    unittest.main()
