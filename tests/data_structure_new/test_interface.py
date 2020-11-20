from coconut.data_structure_new.interface import Interface
from coconut.data_structure_new.model import Model
from coconut.data_structure_new.model_part import ModelPart
import unittest
import numpy as np
import json


class TestInterface(unittest.TestCase):

    def setUp(self):
        model_part_size = 3
        self.parameters = {'interface_a':
                               [{'model_part': 'mp1', 'variables': ['pressure', 'traction']},
                                {'model_part': 'mp2', 'variables': ['density']}]
                           }

        self.model = Model()
        x0 = np.random.rand(3 * model_part_size)
        y0 = np.random.rand(3 * model_part_size)
        z0 = np.random.rand(3 * model_part_size)
        ids = np.arange(0, 3 * model_part_size)
        self.model.create_model_part('mp1', x0[0:model_part_size], y0[0:model_part_size], z0[0:model_part_size],
                                     ids[0:model_part_size])
        self.model.create_model_part('mp2', x0[model_part_size:2 * model_part_size],
                                     y0[model_part_size:2 * model_part_size], z0[model_part_size:2 * model_part_size],
                                     ids[model_part_size: 2 * model_part_size])
        self.model.create_model_part('mp3', x0[2 * model_part_size:3 * model_part_size],
                                     y0[2 * model_part_size:3 * model_part_size],
                                     z0[2 * model_part_size:3 * model_part_size], ids[2*model_part_size: 3 * model_part_size])

        self.interface = Interface(self.parameters['interface_a'], self.model)
        self.scalar_size = 1
        self.vector_size = 3
        self.pressure = np.random.rand(model_part_size, self.scalar_size)
        self.traction = np.random.rand(model_part_size, self.vector_size)
        self.temperature = np.random.rand(model_part_size, self.scalar_size)
        self.density = np.random.rand(model_part_size, self.scalar_size)

    def test_parameters(self):
        parameters = {'mp1': ['pressure']}
        self.assertRaises(TypeError, Interface, parameters, self.model)
        parameters = [{'model_part': 'mp1', 'variables': ['pressure']}, 2]
        self.assertRaises(TypeError, Interface, parameters, self.model)
        parameters = [{'model_part': 'mp1'}]
        self.assertRaises(KeyError, Interface, parameters, self.model)
        parameters = [{'model_part': 'mp1', 'variables': ['pressure'], 'extra': 2}]
        self.assertRaises(KeyError, Interface, parameters, self.model)

        self.assertEqual(self.interface.parameters, self.parameters['interface_a'])

        with self.assertRaises(AttributeError):
            self.interface.parameters = parameters

    def test_model_part_variable_pairs(self):
        ref_result = [('mp1', 'pressure'), ('mp1', 'traction'), ('mp2', 'density')]
        self.assertListEqual(ref_result, self.interface.model_part_variable_pairs)

    def test_set_variable_data(self):
        self.interface.set_variable_data('mp1', 'pressure', self.pressure)
        self.interface.set_variable_data('mp1', 'traction', self.traction)
        self.interface.set_variable_data('mp2', 'density', self.density)
        np.testing.assert_array_equal(self.interface._Interface__data['mp1']['pressure'], self.pressure)
        np.testing.assert_array_equal(self.interface._Interface__data['mp1']['traction'], self.traction)
        np.testing.assert_array_equal(self.interface._Interface__data['mp2']['density'], self.density)






if __name__ == '__main__':
    unittest.main()
