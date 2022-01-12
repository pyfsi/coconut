from coconut import data_structure
from coconut.data_structure.interface import Interface
from coconut.tools import create_instance

import unittest
import numpy as np


class TestPredictorSurrogate(unittest.TestCase):
    
    def setUp(self) -> None:
        m = 10
        dz = 3
        self.a0 = 1
        variable = 'area'
        model_part_name = 'wall'
        self.interface_settings = [{'model_part': model_part_name, 'variables': [variable]}]
    
        # create model and model_part
        model = data_structure.Model()
        ids = np.arange(0, m)
        x0 = np.zeros(m)
        y0 = np.zeros(m)
        z0 = np.arange(0, m * dz, dz)
        model.create_model_part(model_part_name, x0, y0, z0, ids)
    
        a0_array = np.full((m, 1), self.a0)
    
        # create self.interface
        self.interface = Interface(self.interface_settings, model)
        self.interface.set_variable_data(model_part_name, variable, a0_array)
        self.interface_as_array = self.interface.get_interface_data()

    def test_predictor(self):
        s0 = 0.5
        a0 = self.a0
        s1 = 2.5
        p1d = 2.5  # direct prediction
        p1c = 3  # change prediction
        a1 = 2
        s2 = 3.5
        p2d = 3.5  # direct prediction
        p2c = 3  # change prediction
        a2 = 3

        # create predictor directly using surrogate solution
        parameters = {'type': 'predictors.surrogate', 'settings': {'predict_change': False}}
        predictor_surrogate_direct = create_instance(parameters)
        self.interface.set_interface_data(a0 * self.interface_as_array)
        predictor_surrogate_direct.initialize(self.interface)
        self.interface.set_interface_data(s0 * self.interface_as_array)
        predictor_surrogate_direct.update_surrogate(self.interface)

        # create predictor using change in surrogate solution
        parameters = {'type': 'predictors.surrogate'}
        predictor_surrogate_change = create_instance(parameters)
        self.interface.set_interface_data(a0 * self.interface_as_array)
        predictor_surrogate_change.initialize(self.interface)
        self.interface.set_interface_data(s0 * self.interface_as_array)
        predictor_surrogate_change.update_surrogate(self.interface)

        for predictor_surrogate, predictions in zip((predictor_surrogate_direct, predictor_surrogate_change),
                                                    ((p1d, p2d), (p1c, p2c))):
            predictor_surrogate.initialize_solution_step()
            self.interface.set_interface_data(s1 * self.interface_as_array)
            predictor_surrogate.update_surrogate(self.interface)
            self.interface.set_interface_data(a0 * self.interface_as_array)
            prediction = predictor_surrogate.predict(self.interface)
            self.interface.set_interface_data(predictions[0] * self.interface_as_array)
            self.assertEqual(prediction, self.interface)
            self.interface.set_interface_data(a1 * self.interface_as_array)
            predictor_surrogate.update(self.interface)
            predictor_surrogate.finalize_solution_step()

            predictor_surrogate.initialize_solution_step()
            self.interface.set_interface_data(s2 * self.interface_as_array)
            predictor_surrogate.update_surrogate(self.interface)
            self.interface.set_interface_data(a1 * self.interface_as_array)
            prediction = predictor_surrogate.predict(self.interface)
            self.interface.set_interface_data(predictions[1] * self.interface_as_array)
            self.assertEqual(prediction, self.interface)
            self.interface.set_interface_data(a2 * self.interface_as_array)
            predictor_surrogate.update(self.interface)
            predictor_surrogate.finalize_solution_step()

    def test_errors(self):
        # create predictor
        parameters = {'type': 'predictors.surrogate'}
        predictor_surrogate = create_instance(parameters)
        predictor_surrogate.initialize(self.interface)

        # error if no update
        with self.assertRaises(Exception):
            predictor_surrogate.initialize_solution_step()
            predictor_surrogate.finalize_solution_step()

        # error if no surrogate update
        with self.assertRaises(Exception):
            predictor_surrogate.initialize_solution_step()
            predictor_surrogate.update(self.interface)
            predictor_surrogate.finalize_solution_step()

        # error if updated twice
        with self.assertRaises(Exception):
            predictor_surrogate.initialize_solution_step()
            _ = predictor_surrogate.update(self.interface)
            _ = predictor_surrogate.update(self.interface)

        # error if surrogate updated twice
        with self.assertRaises(Exception):
            predictor_surrogate.initialize_solution_step()
            _ = predictor_surrogate.update_surrogate(self.interface)
            _ = predictor_surrogate.update_surrogate(self.interface)

        # error if prediction before surrogate update
        with self.assertRaises(Exception):
            predictor_surrogate.initialize_solution_step()
            _ = predictor_surrogate.predict(self.interface)

        # error if prediction after update
        with self.assertRaises(Exception):
            predictor_surrogate.initialize_solution_step()
            _ = predictor_surrogate.update_surrogate(self.interface)
            _ = predictor_surrogate.update(self.interface)
            _ = predictor_surrogate.predict(self.interface)
            predictor_surrogate.finalize_solution_step()


if __name__ == '__main__':
    unittest.main()
