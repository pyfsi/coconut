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
                            {'model_part': 'mp2', 'variables': ['density']},
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
        self.temperature = np.random.rand(model_part_size, self.scalar_size)
        self.density = np.random.rand(model_part_size, self.scalar_size)
        self.interface_data = np.random.rand(model_part_size * 5)

    def test_instantiation(self):
        parameters = {'mp1': ['pressure']}
        self.assertRaises(TypeError, Interface, parameters, self.model)
        parameters = [{'model_part': 'mp1', 'variables': ['pressure']}, 2]
        self.assertRaises(TypeError, Interface, parameters, self.model)
        parameters = [{'model_part': 'mp1'}]
        self.assertRaises(KeyError, Interface, parameters, self.model)
        parameters = [{'model_part': 'mp1', 'variables': ['pressure'], 'extra': 2}]
        self.assertRaises(KeyError, Interface, parameters, self.model)
        parameters = [{'model_part': 'mp1', 'variables': 'pressure'}]
        self.assertRaises(TypeError, Interface, parameters, self.model)

        self.assertEqual(self.interface.parameters, self.parameters['interface_a'])

    def test_model_part_variable_pairs(self):
        ref_result = [('mp1', 'pressure'), ('mp1', 'traction'), ('mp2', 'density')]
        self.assertListEqual(ref_result, self.interface.model_part_variable_pairs)

    def test_size(self):
        self.assertEqual(self.interface.size, 5 * self.model_part_size)

    def test_attribute_change(self):
        with self.assertRaises(AttributeError):
            self.interface.model_part_variable_pairs = ()
        with self.assertRaises(AttributeError):
            parameters = {'mp1': ['pressure']}
            self.interface.parameters = parameters
        with self.assertRaises(AttributeError):
            self.interface.size = 0

    def test_copy(self):
        self.interface.set_interface_data(self.interface_data)
        interface_copy = self.interface.copy()
        self.assertIsNot(self.interface, interface_copy)
        # interface_copy has the same values
        self.assertEqual(self.interface.model_part_variable_pairs, interface_copy.model_part_variable_pairs)
        np.testing.assert_array_equal(self.interface.get_interface_data(), interface_copy.get_interface_data())
        # interface_copy not affected by change in self.interface
        interface_data = self.interface.get_interface_data()
        self.interface.set_interface_data(interface_data * 2)
        np.testing.assert_array_equal(interface_copy.get_interface_data(), interface_data)

    def test_get_variable_data(self):
        self.assertRaises(KeyError, self.interface.get_variable_data, 'mp2', 'pressure')

        model_part_variable_pair = ('mp1', 'pressure')
        variable_data = self.pressure
        # returns copy
        self.interface.set_variable_data(*model_part_variable_pair, variable_data)
        variable_data = self.interface.get_variable_data(*model_part_variable_pair)
        variable_data *= 2
        np.testing.assert_array_equal(self.interface.get_variable_data(*model_part_variable_pair), variable_data / 2)

    def test_set_variable_data(self):
        self.interface.set_variable_data('mp1', 'pressure', self.pressure)
        self.interface.set_variable_data('mp1', 'traction', self.traction)
        self.interface.set_variable_data('mp2', 'density', self.density)
        np.testing.assert_array_equal(self.interface._Interface__data['mp1']['pressure'], self.pressure)
        np.testing.assert_array_equal(self.interface._Interface__data['mp1']['traction'], self.traction)
        np.testing.assert_array_equal(self.interface._Interface__data['mp2']['density'], self.density)
        # input is array with correct shape
        self.assertRaises(ValueError, self.interface.set_variable_data, 'mp1', 'pressure', self.traction)
        self.assertRaises(ValueError, self.interface.set_variable_data, 'mp1', 'pressure', self.pressure.flatten())
        self.assertRaises(ValueError, self.interface.set_variable_data, 'mp1', 'pressure', list(self.pressure))
        # copy is made of input
        self.interface.set_variable_data('mp1', 'pressure', self.pressure)
        self.pressure *= 2
        np.testing.assert_array_equal(self.interface._Interface__data['mp1']['pressure'], self.pressure / 2)

    def test_get_interface_data(self):
        # output has correct size
        self.assertEqual(self.interface.get_interface_data().size, self.interface.size)
        # correct output from variable data
        self.interface.set_variable_data('mp1', 'pressure', self.pressure)
        self.interface.set_variable_data('mp1', 'traction', self.traction)
        self.interface.set_variable_data('mp2', 'density', self.density)
        interface_data = np.concatenate((self.pressure.flatten(), self.traction.flatten(), self.density.flatten()))
        np.testing.assert_equal(self.interface.get_interface_data(), interface_data)
        # correct output from interface data
        self.interface.set_interface_data(self.interface_data)
        np.testing.assert_equal(self.interface.get_interface_data(), self.interface_data)
        # returns copy
        self.interface.get_interface_data() * 2
        np.testing.assert_equal(self.interface.get_interface_data(), self.interface_data)

    def test_set_interface_data(self):
        self.interface.set_interface_data(self.interface_data)
        np.testing.assert_array_equal(self.interface._Interface__data['mp1']['pressure'].flatten(),
                                      self.interface_data[:self.scalar_size * self.model_part_size])
        np.testing.assert_array_equal(self.interface._Interface__data['mp1']['traction'].flatten(),
                                      self.interface_data[self.scalar_size * self.model_part_size:
                                                          (self.scalar_size + self.vector_size) * self.model_part_size])
        np.testing.assert_array_equal(self.interface._Interface__data['mp2']['density'].flatten(),
                                      self.interface_data[(self.scalar_size + self.vector_size)
                                                          * self.model_part_size:])
        # input is array with correct shape
        self.assertRaises(ValueError, self.interface.set_interface_data, self.pressure)
        self.assertRaises(ValueError, self.interface.set_interface_data, self.interface_data.reshape(-1, 1))
        self.assertRaises(ValueError, self.interface.set_interface_data, list(self.interface_data))
        # copy is made of input
        self.interface.set_interface_data(self.interface_data)
        self.interface_data *= 2
        np.testing.assert_array_equal(self.interface.get_interface_data(), self.interface_data / 2)

    def test_norm(self):
        self.interface.set_interface_data(self.interface_data)
        norm = np.linalg.norm(self.interface_data)
        self.assertEqual(self.interface.norm(), norm)

    def create_test_interfaces(self):
        interface_data1 = np.random.rand(self.interface.size)
        interface_data2 = np.random.rand(self.interface.size)
        interface1 = self.interface.copy()
        interface1.set_interface_data(interface_data1)
        interface2 = self.interface.copy()
        interface2.set_interface_data(interface_data2)
        return interface_data1, interface_data2, interface1, interface2

    def test_add(self):
        interface_data1, interface_data2, interface1, interface2 = self.create_test_interfaces()

        interface_sum = interface1 + interface2
        self.interface.set_interface_data(interface_data1 + interface_data2)
        self.assertEqual(interface_sum, self.interface)
        interface3 = interface1.copy()
        interface3 += interface2
        self.assertEqual(interface3, self.interface)

        for number in (np.random.rand(1), float(np.random.rand(1)), int(10 * np.random.rand(1))):
            interface_sum = interface1 + number
            self.interface.set_interface_data(interface_data1 + number)
            self.assertEqual(interface_sum, self.interface)
            interface_sum = number + interface1
            self.assertEqual(interface_sum, self.interface)
            interface3 = interface1.copy()
            interface3 += number
            self.assertEqual(interface3, self.interface)

        for other in ('a', True):
            with self.assertRaises(TypeError):
                _ = interface1 + other
            with self.assertRaises(TypeError):
                _ = other + interface1
            with self.assertRaises(TypeError):
                interface1 += other

    def test_sub(self):
        interface_data1, interface_data2, interface1, interface2 = self.create_test_interfaces()

        interface_sum = interface1 - interface2
        self.interface.set_interface_data(interface_data1 - interface_data2)
        self.assertEqual(interface_sum, self.interface)
        interface3 = interface1.copy()
        interface3 -= interface2
        self.assertEqual(interface3, self.interface)

        for number in (np.random.rand(1), float(np.random.rand(1)), int(10 * np.random.rand(1))):
            interface_sum = interface1 - number
            self.interface.set_interface_data(interface_data1 - number)
            self.assertEqual(interface_sum, self.interface)
            interface_sum = number - interface1
            self.assertEqual(interface_sum, self.interface)
            interface3 = interface1.copy()
            interface3 -= number
            self.assertEqual(interface3, self.interface)

        for other in ('a', True):
            with self.assertRaises(TypeError):
                _ = interface1 - other
            with self.assertRaises(TypeError):
                _ = other - interface1
            with self.assertRaises(TypeError):
                interface1 -= other

    def test_mul(self):
        interface_data1, interface_data2, interface1, interface2 = self.create_test_interfaces()

        for number in (np.random.rand(1), float(np.random.rand(1)), int(10 * np.random.rand(1))):
            interface_sum = interface1 * number
            self.interface.set_interface_data(interface_data1 * number)
            self.assertEqual(interface_sum, self.interface)
            interface_sum = number * interface1
            self.assertEqual(interface_sum, self.interface)
            interface3 = interface1.copy()
            interface3 *= number
            self.assertEqual(interface3, self.interface)

        for other in ('a', True, interface2):
            with self.assertRaises(TypeError):
                _ = interface1 * other
            with self.assertRaises(TypeError):
                _ = other * interface1
            with self.assertRaises(TypeError):
                interface1 *= other

    def test_truediv(self):
        interface_data1, interface_data2, interface1, interface2 = self.create_test_interfaces()

        for number in (np.random.rand(1), float(np.random.rand(1)), int(10 * np.random.rand(1))):
            interface_sum = interface1 / number
            self.interface.set_interface_data(interface_data1 / number)
            self.assertEqual(interface_sum, self.interface)
            with self.assertRaises(TypeError):
                _ = number / interface1
            interface3 = interface1.copy()
            interface3 /= number
            self.assertEqual(interface3, self.interface)

        for other in ('a', True, interface2):
            with self.assertRaises(TypeError):
                _ = interface1 / other
            with self.assertRaises(TypeError):
                _ = other / interface1
            with self.assertRaises(TypeError):
                interface1 /= other


if __name__ == '__main__':
    unittest.main()
