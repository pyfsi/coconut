from coconut.data_structure.interface import Interface
from coconut.data_structure.model import Model

import unittest
import numpy as np


class TestInterface(unittest.TestCase):

    def setUp(self):
        self.model_part_size = model_part_size = 3
        self.parameters = {'interface_a':
                           [
                            {'model_part': 'mp1', 'variables': ['pressure', 'traction']},
                            {'model_part': 'mp2', 'variables': ['displacement']},
                            ]
                           }

        self.model = Model()
        x0 = np.random.rand(3 * model_part_size)
        y0 = np.random.rand(3 * model_part_size)
        z0 = np.random.rand(3 * model_part_size)
        ids = np.arange(0, model_part_size)
        self.model.create_model_part('mp1', x0[0:model_part_size], y0[0:model_part_size], z0[0:model_part_size], ids)
        self.model.create_model_part('mp2', x0[model_part_size:2 * model_part_size],
                                     y0[model_part_size:2 * model_part_size],
                                     z0[model_part_size:2 * model_part_size], ids)
        self.model.create_model_part('mp3', x0[2 * model_part_size:3 * model_part_size],
                                     y0[2 * model_part_size:3 * model_part_size],
                                     z0[2 * model_part_size:3 * model_part_size], ids)

        self.interface = Interface(self.parameters['interface_a'], self.model)
        self.scalar_size = 1
        self.vector_size = 3
        self.pressure = np.random.rand(model_part_size, self.scalar_size)
        self.traction = np.random.rand(model_part_size, self.vector_size)
        self.displacement = np.random.rand(model_part_size, self.vector_size)
        self.displacement[:,2] = 0
        self.density = np.random.rand(model_part_size, self.scalar_size)
        self.interface_data = np.random.rand(model_part_size * 7)

    def test_model_part_variable_pairs(self):
        ref_result = [('mp1', 'pressure'), ('mp1', 'traction'), ('mp2', 'displacement')]
        self.assertListEqual(ref_result, self.interface.model_part_variable_pairs)

    def test_get_variable_data(self):
        self.assertRaises(KeyError, self.interface.get_variable_data, 'mp2', 'pressure')

        model_part_variable_pair = ('mp2', 'displacement')
        variable_data = self.displacement
        #
        self.interface.set_variable_data(*model_part_variable_pair, variable_data)
        variable_data = self.interface.get_variable_data(*model_part_variable_pair)

        print(variable_data)
        print(variable_data[0])

        r = self.model.get_model_part('mp1')
        x_out = r.x0
        y_out =r.y0
        print(x_out)
        print(y_out)
        variable_data[0]= np.add(variable_data[0], x_out)
        variable_data[1]= np.add(variable_data[1], y_out)

        print("mp_from")
        print(variable_data)
        ids = np.arange(0,3)
        self.model.create_model_part('new_mp', variable_data[:,0], variable_data[:,1], variable_data[:,2], ids)
        l = self.model.get_model_part('new_mp')
        x = l.x0
        y = l.y0
        z = l.z0
        print(x,y,z)

    def test_set_variable_data(self):
        self.interface.set_variable_data('mp1', 'pressure', self.pressure)
        self.interface.set_variable_data('mp1', 'traction', self.traction)
        self.interface.set_variable_data('mp2', 'displacement', self.displacement)
        np.testing.assert_array_equal(self.interface._Interface__data['mp1']['pressure'], self.pressure)
        np.testing.assert_array_equal(self.interface._Interface__data['mp1']['traction'], self.traction)
        np.testing.assert_array_equal(self.interface._Interface__data['mp2']['displacement'], self.displacement)
        # input is array with correct shape
        self.assertRaises(ValueError, self.interface.set_variable_data, 'mp1', 'pressure', self.traction)
        self.assertRaises(ValueError, self.interface.set_variable_data, 'mp1', 'pressure', self.pressure.flatten())
        self.assertRaises(ValueError, self.interface.set_variable_data, 'mp1', 'pressure', list(self.pressure))
        # copy is made of input
        self.interface.set_variable_data('mp1', 'pressure', self.pressure)
        self.pressure *= 2
        np.testing.assert_array_equal(self.interface._Interface__data['mp1']['pressure'], self.pressure / 2)





if __name__ == '__main__':
    unittest.main()
