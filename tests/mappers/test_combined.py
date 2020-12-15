from coconut import data_structure
import unittest
from coconut.coupling_components.tools import create_instance

import numpy as np
import os

class TestMapperCombined(unittest.TestCase):
    def test_mapper_combined(self):
        parameter_file_name = os.path.join(os.path.dirname(__file__), 'test_combined.json')
        with open(parameter_file_name, 'r') as parameter_file:
            parameters = data_structure.Parameters(parameter_file.read())

        # compare 3 mappers: nearest, combined_a, combined_b
        if True:
            # create model_part_from
            var_from = vars(data_structure)["TEMPERATURE"]
            model_from = data_structure.Model()
            model_part_from = model_from.CreateModelPart('wall_from')
            model_part_from.AddNodalSolutionStepVariable(var_from)

            for i in range(100):
                node = model_part_from.CreateNewNode(i, i, i ** 2, i ** 3)
                node.SetSolutionStepValue(var_from, 0, np.random.rand())

            # create model_part_to
            var_to = vars(data_structure)["PRESSURE"]
            model_to = data_structure.Model()
            model_part_to = model_to.CreateModelPart('wall_to')
            model_part_to.AddNodalSolutionStepVariable(var_to)

            for i in range(100):
                model_part_to.CreateNewNode(i, (99 - i), (99 - i) ** 2, (99 - i) ** 3)

            # create mappers, get output data
            data = []
            for mapper_name in ['mapper_nearest', 'mapper_combined_a', 'mapper_combined_b']:
                mapper = create_instance(parameters[mapper_name])
                mapper.initialize(model_part_from, model_part_to)
                mapper((model_part_from, var_from), (model_part_to, var_to))

                values = []
                for i_to, node in enumerate(model_part_to.Nodes):
                    values.append(node.GetSolutionStepValue(var_to))
                data.append(values)

            # check output data
            self.assertListEqual(data[0], data[1])
            self.assertListEqual(data[0], data[2])


if __name__ == '__main__':
    unittest.main()
