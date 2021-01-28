from coconut.data_structure.interface import Interface
from coconut.data_structure.model import Model

import unittest
import numpy as np


class TestInterfaceMapper(unittest.TestCase):

    def __init__(self,n_from,n_to):

        self.parameters = {'type': 'mappers.interface_new_mp',
                           'settings':
                               {
                                'type': 'mappers.linear',
                                'settings':
                                    {'directions':
                                        [
                                        "X"
                                        ]
                                    }
                                }
                           }

        model_part_size = 3
        self.model = Model()
        self.settings = self.parameters['settings']
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




if __name__ == '__main__':
    unittest.main()
