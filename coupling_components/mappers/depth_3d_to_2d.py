from coconut.coupling_components.mappers.depth_2d_to_3d import MapperDepth2DTo3D
from coconut.coupling_components.mappers.transformer import MapperTransformer
from coconut.data_structure import variables_dimensions

import numpy as np


def create(parameters):
    return MapperDepth3DTo2D(parameters)


class MapperDepth3DTo2D(MapperDepth2DTo3D):
    def initialize(self, model, model_part_name_in, model_part_name_out, forward):
        if not forward:
            super().initialize(model, model_part_name_in, model_part_name_out, not forward)
        else:
            raise NotImplementedError('Forward Initialization not implemented for MapperDepth3DTo2D.')
        # switch from and to variables
        n_from = self.n_to
        n_to = self.n_from
        self.n_from = n_from
        self.n_to = n_to

    def __call__(self, args_from, args_to):
        MapperTransformer.__call__(self, args_from, args_to)

        interface_from, mp_name_from, var = args_from
        interface_to, mp_name_to, _ = args_to

        dimensions = variables_dimensions[var]
        data_from = interface_from.get_variable_data(mp_name_from, var)
        data_to = np.zeros((self.n_to, dimensions))
        for i_d in range(self.n_depth):
            i_start = i_d * self.n_to
            i_end = (i_d + 1) * self.n_to
            data_to += data_from[i_start:i_end] / self.n_depth
        if dimensions == 3:
            data_to[:, self.dir_depth] *= 0
        elif dimensions != 1:
            raise NotImplementedError(
                f'MapperDepth2DTo3D not implemented for variable of dimension {dimensions}')
        interface_to.set_variable_data(mp_name_to, var, data_to)
