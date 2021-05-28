from coconut.coupling_components.mappers.transformer import MapperTransformer
from coconut.coupling_components.mappers.axisymmetric_2d_to_3d_mod import MapperAxisymmetric2DTo3DMod
from coconut.data_structure import variables_dimensions

import numpy as np


def create(parameters):
    return MapperAxisymmetric3DTo2DMod(parameters)


class MapperAxisymmetric3DTo2DMod(MapperAxisymmetric2DTo3DMod):
    def initialize(self, model, model_part_name_in, model_part_name_out, forward):
        if not forward:
            super().initialize(model, model_part_name_in, model_part_name_out, not(forward))

            n_from = self.n_to
            n_to = self.n_from
            self.n_from = n_from
            self.n_to = n_to

        else:
            raise NotImplementedError('Forward Initialization not implemented for MapperAxisymmetric3DTo2D.')

    def __call__(self, args_from, args_to):
        MapperTransformer.__call__(self, args_from, args_to)

        interface_from, mp_name_from, var = args_from
        interface_to, mp_name_to, _ = args_to

        dimensions = variables_dimensions[var]
        data_from = interface_from.get_variable_data(mp_name_from, var)
        data_to = np.zeros((self.n_to, dimensions))
        if dimensions == 1:
            for i_t in range(self.n_t):
                i_start = i_t * self.n_to
                i_end = (i_t + 1) * self.n_to
                data_to += data_from[i_start:i_end] / self.n_t
        elif dimensions == 3:

            r = (np.cos(self.theta) * data_from[:, self.dir_r] +
                 np.sin(self.theta) * data_from[:, self.dir_3d])
            data_from[:, self.dir_r] = r
            data_from[:, self.dir_3d] *= 0
            for i_t in range(self.n_t):
                i_start = i_t * self.n_to
                i_end = (i_t + 1) * self.n_to
                data_to += data_from[i_start:i_end] / self.n_t
        else:
            raise NotImplementedError(
                f'MapperAxisymmetric2DTo3D not implemented for variable of dimension {dimensions}')
        interface_to.set_variable_data(mp_name_to, var, data_to)
