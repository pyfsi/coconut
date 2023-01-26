from coconut import data_structure
from coconut.data_structure.interface import Interface
from coconut.tools import create_instance

import unittest
import numpy as np
import copy


class TestPredictorSurrogate(unittest.TestCase):
    def setUp(self):
        m = 10
        dz = 3
        self.a0 = 1
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

        a0_array = np.full((m, 1), self.a0)

        # create interface
        self.interface = Interface(interface_settings, model)
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

    def test_change_settings_on_restart(self):
        # test if changing settings upon restart works correctly

        s0 = 0.5
        s1 = 2.5
        a1 = 2
        s2 = 3.5
        a2 = 3
        s3 = 0
        a3 = -0.5

        # create predictor directly using surrogate solution
        parameters = {'type': 'predictors.surrogate', 'settings': {'predict_change': False}}
        predictor = create_instance(parameters)

        predictor.initialize(self.interface)
        predictor.update_surrogate(s0 * self.interface)

        def do_time_step(s, a):
            predictor.initialize_solution_step()
            predictor.update_surrogate(s * self.interface)
            predictor.update(a * self.interface)
            predictor.finalize_solution_step()

        do_time_step(s1, a1)

        self.restart_data = copy.deepcopy(predictor.save_restart_data())
        self.parameters_old = parameters.copy()
        dataprev = predictor.dataprev.copy()
        surrogate_data = predictor.surrogate_data.copy()
        predictor.finalize()

        # create restarted predictor
        predict_change_new = True
        parameters = {'type': 'predictors.surrogate', 'settings': {'predict_change': predict_change_new}}
        predictor = create_instance(parameters)
        predictor.check_restart_data(self.parameters_old)
        predictor.initialize(self.interface)
        predictor.restart(self.restart_data)

        self.assertEqual(predictor.predict_change, predict_change_new)

        np.testing.assert_allclose(predictor.dataprev, dataprev)
        np.testing.assert_allclose(predictor.surrogate_data, surrogate_data)

        predictor.update_surrogate(s2 * self.interface)
        self.assertEqual(predictor.predict(self.interface), 3 * self.interface)
        predictor.update(a2 * self.interface)

        do_time_step(s3, a3)

        predictor.finalize()


if __name__ == '__main__':
    unittest.main()
