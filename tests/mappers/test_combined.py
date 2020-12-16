from coconut import data_structure
from coconut.coupling_components.tools import create_instance

import unittest
import numpy as np
import json

class TestMapperCombined(unittest.TestCase):
    def test_mapper_combined(self):
        with open('mappers/test_combined.json') as parameter_file:
            parameters = json.load(parameter_file)

        # compare 3 mappers: nearest, combined_a, combined_b
        for var in ['pressure', 'displacement']:
            var = 'pressure'
            mp_name_from = 'wall_from'
            mp_name_to = 'wall_to'
            model = data_structure.Model()

            n = 100
            tmp = np.linspace(0, 1, n)
            x, y, z = tmp, tmp ** 1.1, tmp ** 1.2
            v_from = np.random.rand(n, 1)
            mp_from = model.create_model_part(mp_name_from, x, y, z, np.arange(n))
            mp_to = model.create_model_part(mp_name_to, np.flip(x), np.flip(y), np.flip(z), np.arange(n))

            parameters_from = [{'model_part': mp_name_from, 'variables': [var]}]
            int_from = data_structure.Interface(parameters_from, model)
            int_from.set_variable_data(mp_name_from, var, v_from)
            parameters_to = [{'model_part': mp_name_to, 'variables': [var]}]
            int_to = data_structure.Interface(parameters_to, model)

            # create mappers, get output data
            data = []
            for mapper_name in ['mapper_nearest', 'mapper_combined_a', 'mapper_combined_b']:
                mapper = create_instance(parameters[mapper_name])
                mapper.initialize(mp_from, mp_to)
                mapper((int_from, mp_name_from, var), (int_to, mp_name_to, var))
                data.append(int_to.get_variable_data(mp_name_to, var))

            # check output data
            np.testing.assert_array_equal(data[0], data[1])
            np.testing.assert_array_equal(data[0], data[2])


if __name__ == '__main__':
    unittest.main()
