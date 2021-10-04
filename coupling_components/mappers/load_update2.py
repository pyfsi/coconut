from coconut.coupling_components.mappers.transformer import MapperTransformer
from coconut.data_structure import variables_dimensions
from scipy.spatial import cKDTree
from coconut import data_structure
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt

import numpy as np

def create(parameters):
    return Mapper_load_update2(parameters)


class Mapper_load_update2(MapperTransformer):
    def __init__(self, parameters):
        super().__init__(parameters)

        self.v_min = self.settings['coords_min']
        if type(self.v_min) != float:
            raise TypeError('coords_min must be a float')

        self.v_max = self.settings['coords_max']
        if type(self.v_max) != float:
            raise TypeError('coords_max must be a float')

    def initialize(self, model, model_part_name_in, model_part_name_out, forward):
        if not forward:
            super().initialize()

            self.mp_input_to = model.get_model_part(model_part_name_in)

            n_in = self.mp_input_to.size
            n_out = n_in
            self.n_from = n_out
            self.n_to = n_in

            self.coords_in = np.column_stack((self.mp_input_to.x0, self.mp_input_to.y0, self.mp_input_to.z0))
            ids_out = np.arange(n_out)


            self.mp_out = model.create_model_part(model_part_name_out, self.coords_in[:, 0],
                                    self.coords_in[:, 1], self.coords_in[:, 2], ids_out)

        else:
            raise NotImplementedError('Backward Initialization not implemented for LoadUpdate')

    def __call__(self, args_from, args_to):
        super().__call__(args_from, args_to)

        interface_from, mp_name_from, var = args_from
        interface_to, mp_name_to, _ = args_to

        self.v_min = self.settings['coords_min']
        self.v_max = self.settings['coords_max']

        dimensions = variables_dimensions[var]
        self.data_from = interface_from.get_variable_data(mp_name_from, var)
        coord_from = interface_from.get_model_part(mp_name_from)
        # print("coord.x0 mapper load")
        # print(coord_from.y0)
        # print("data_from")
        # print(self.data_from)
        # y = self.data_from[:, 1]
        # print(y)
        #
        if dimensions == 1:
            data_to = np.zeros((self.mp_input_to.size, 1))
            for i in range(self.mp_input_to.size):
                if self.mp_input_to.y0[i] < self.v_min or self.mp_input_to.y0[i] > self.v_max:
                    data_to[i] = 0
                else:
                     data_to[i] = self.data_from[i]
            y = data_to
            # print("coord.x0 mapper load")
            # print(coord_from.y0)
            # print(y)
            # plt.scatter(coord_from.y0, y, color='b', label="input_to")
            # plt.show()
        elif dimensions == 3:

            data_to = np.zeros((self.mp_input_to.size, 3))
            for j in range(self. mp_input_to.size):
                if self.mp_input_to.x0[j] < self.v_min or self.mp_input_to.x0[j] > self.v_max:
                    data_to[j,0] = 0
                    data_to[j,1] = 0
                    data_to[j,2] = 0
                else:
                    data_to[j,0] = self.data_from[j,0]
                    data_to[j, 1] = self.data_from[j,1]
                    data_to[j, 2] = self.data_from[j,2]

        else:
            raise NotImplementedError(
                f'MapperUpdateLoad not implemented for variable of dimension {dimensions}')


        interface_to.set_variable_data(mp_name_to, var, data_to)