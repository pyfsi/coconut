from coconut import data_structure
from coconut.data_structure.interface import Interface
from coconut.tools import create_instance

import unittest
import numpy as np
import copy


class TestPredictor(unittest.TestCase):
    def setUp(self):
        self.m = 10
        dz = 3
        a0 = 1
        variable = 'area'
        model_part_name = 'wall'
        interface_settings = [{'model_part': model_part_name, 'variables': [variable]}]

        # create model and model_part
        model = data_structure.Model()
        ids = np.arange(0, self.m)
        x0 = np.zeros(self.m)
        y0 = np.zeros(self.m)
        z0 = np.arange(0, self.m * dz, dz)
        model.create_model_part(model_part_name, x0, y0, z0, ids)

        a0_array = np.full((self.m, 1), a0)

        # create interface
        self.interface = Interface(interface_settings, model)
        self.interface.set_variable_data(model_part_name, variable, a0_array)

        # create predictor
        self.parameters = {'type': 'predictors.cubic'}
        self.predictor = create_instance(self.parameters)

    def test_predictor(self):
        a1 = 2
        a2 = 3
        a3 = 4
        a4 = 5

        interface = self.interface
        interface_as_array = interface.get_interface_data()

        self.predictor.initialize(interface)

        def do_time_step(a):
            self.predictor.initialize_solution_step()
            interface.set_interface_data(a * interface_as_array)
            self.predictor.update(interface)
            self.predictor.finalize_solution_step()

        for a in (a1, a2, a3):
            do_time_step(a)

        self.predictor.initialize_solution_step()
        prediction_linear = self.predictor.linear(interface).get_interface_data()
        prediction_quadratic = self.predictor.quadratic(interface).get_interface_data()
        prediction_cubic = self.predictor.cubic(interface).get_interface_data()
        for i in range(self.m):
            self.assertAlmostEqual(a4, prediction_linear[i])
            self.assertAlmostEqual(a4, prediction_quadratic[i])
            self.assertAlmostEqual(a4, prediction_cubic[i])

        # error if no update
        with self.assertRaises(Exception):
            self.predictor.initialize_solution_step()
            self.predictor.finalize_solution_step()

        # error if updated twice
        with self.assertRaises(Exception):
            self.predictor.initialize_solution_step()
            _ = self.predictor.update(interface)
            _ = self.predictor.update(interface)

        # error if prediction after update
        with self.assertRaises(Exception):
            self.predictor.initialize_solution_step()
            _ = self.predictor.update(interface)
            _ = self.predictor.predict(interface)

    def test_change_settings_on_restart(self):
        # test if changing settings upon restart works correctly

        a1 = 2
        a2 = 3
        a3 = 4
        a4 = 5
        a5 = 2

        interface = self.interface
        interface_as_array = interface.get_interface_data()

        self.predictor.initialize(interface)

        def do_time_step(a):
            self.predictor.initialize_solution_step()
            interface.set_interface_data(a * interface_as_array)
            self.predictor.update(interface)
            self.predictor.finalize_solution_step()

        # run for 2 time steps
        for a in (a1, a2, a3):
            do_time_step(a)

        self.restart_data = copy.deepcopy(self.predictor.save_restart_data())
        self.parameters_old = self.parameters.copy()
        dataprev = self.predictor.dataprev.copy()
        self.predictor.finalize()

        # create restarted predictor
        self.parameters = {'type': 'predictors.linear'}
        self.predictor = create_instance(self.parameters)
        self.predictor.check_restart_data(self.parameters_old)
        self.predictor.initialize(interface)
        self.predictor.restart(self.restart_data)

        np.testing.assert_allclose(self.predictor.dataprev, [dataprev[0], dataprev[1]])
        np.testing.assert_allclose(self.predictor.predict(self.interface).get_interface_data(), 5 * interface_as_array)

        # run for 2 time steps
        for a in (a4, a5):
            do_time_step(a)

        self.predictor.finalize()


if __name__ == '__main__':
    unittest.main()
