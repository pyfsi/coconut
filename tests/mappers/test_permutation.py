from coconut import data_structure
from coconut.coupling_components.tools import create_instance
from coconut.coupling_components.mappers.permutation import MapperPermutation

import unittest
import numpy as np


class TestMapperPermutation(unittest.TestCase):

    def setUp(self):
        self.parameters = {'type': 'mappers.permutation',
                           'settings': {'permutation': [2, 0, 1]}}

    def test_instantiation(self):
        self.parameters['settings']['permutation'] = "some string"
        self.assertRaises(TypeError, MapperPermutation, self.parameters)

        self.parameters['settings']['permutation'] = [0, 1]
        self.assertRaises(ValueError, MapperPermutation, self.parameters)

        self.parameters['settings']['permutation'] = [1, 2, 3]
        self.assertRaises(ValueError, MapperPermutation, self.parameters)

        self.parameters['settings']['permutation'] = [0, 1, 1]
        self.assertRaises(ValueError, MapperPermutation, self.parameters)

    def test_initialize(self):
        mp_name_in = 'wall_in'
        mp_name_out_f = 'wall_out_f'
        mp_name_out_b = 'wall_out_b'

        n = 10
        x, y, z = np.random.rand(n), np.random.rand(n), np.random.rand(n)

        model = data_structure.Model()
        model.create_model_part(mp_name_in, x, y, z, np.arange(n))

        # model_part_from given
        mapper = create_instance(self.parameters)
        mapper.initialize(model, mp_name_in, mp_name_out_f, forward=True)
        mp_out = model.get_model_part(mp_name_out_f)
        coords_out = np.column_stack((mp_out.x0, mp_out.y0, mp_out.z0))
        np.testing.assert_array_equal(coords_out, np.column_stack((z, x, y)))

        # model_part_to given
        mapper = create_instance(self.parameters)
        mapper.initialize(model, mp_name_in, mp_name_out_b, forward=False)
        mp_out = model.get_model_part(mp_name_out_b)
        coords_out = np.column_stack((mp_out.x0, mp_out.y0, mp_out.z0))
        np.testing.assert_array_equal(coords_out, np.column_stack((y, z, x)))

    def test_call_1d_var(self):
        var = 'pressure'
        mp_name_in = 'wall_in'
        mp_name_out = 'wall_out'

        n = 10
        x, y, z = np.random.rand(n), np.random.rand(n), np.random.rand(n)
        v_in = np.random.rand(n, 1)

        model = data_structure.Model()
        model.create_model_part(mp_name_in, x, y, z, np.arange(n))
        parameters_in= [{'model_part': mp_name_in, 'variables': [var]}]
        interface_in = data_structure.Interface(parameters_in, model)
        interface_in.set_variable_data(mp_name_in, var, v_in)

        mapper = create_instance(self.parameters)
        mapper.initialize(model, mp_name_in, mp_name_out, forward=True)

        parameters_out = [{'model_part': mp_name_out, 'variables': [var]}]
        interface_out = data_structure.Interface(parameters_out, model)
        mapper((interface_in, mp_name_in, var),
               (interface_out, mp_name_out, var))

        v_out = interface_out.get_variable_data(mp_name_out, var)
        np.testing.assert_array_equal(v_in, v_out)

    def test_call_3d_var(self):
        var = 'displacement'
        mp_name_in = 'wall_in'
        mp_name_out = 'wall_out'

        n = 10
        x, y, z = np.random.rand(n), np.random.rand(n), np.random.rand(n)
        v_in = np.random.rand(n, 3)

        model = data_structure.Model()
        model.create_model_part(mp_name_in, x, y, z, np.arange(n))
        parameters_in= [{'model_part': mp_name_in, 'variables': [var]}]
        interface_in = data_structure.Interface(parameters_in, model)
        interface_in.set_variable_data(mp_name_in, var, v_in)

        mapper = create_instance(self.parameters)
        mapper.initialize(model, mp_name_in, mp_name_out, forward=True)

        parameters_out = [{'model_part': mp_name_out, 'variables': [var]}]
        interface_out = data_structure.Interface(parameters_out, model)
        mapper((interface_in, mp_name_in, var),
               (interface_out, mp_name_out, var))

        v_out = interface_out.get_variable_data(mp_name_out, var)
        np.testing.assert_array_equal(v_in[:, 0], v_out[:, 1])
        np.testing.assert_array_equal(v_in[:, 1], v_out[:, 2])
        np.testing.assert_array_equal(v_in[:, 2], v_out[:, 0])


if __name__ == '__main__':
    unittest.main()
