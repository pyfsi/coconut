from coconut.coupling_components.mappers.transformer import MapperTransformer
from coconut.data_structure import variables_dimensions

import numpy as np


def create(parameters):
    return MapperDepth2DTo3D(parameters)


class MapperDepth2DTo3D(MapperTransformer):
    def __init__(self, parameters):
        super().__init__(parameters)

        # get depth direction
        dirs = ['x', 'y', 'z']
        if self.settings['direction_depth'] not in dirs:
            raise ValueError(f'Invalid depth direction {self.settings["direction_depth"]}')
        self.dir_depth = dirs.index(self.settings['direction_depth'])
        self.dir_other = [0, 1, 2]
        self.dir_other.pop(self.dir_depth)

        # get coordinates in depth direction
        self.depth_coords = self.settings['coordinates_depth']
        if type(self.depth_coords) != list:
            raise ValueError(f'"coordinates_depth" must be a list')

        # get number of nodes in depth direction
        self.n_depth = len(self.depth_coords)

        self.n_from = self.n_to = None

    def initialize(self, model, model_part_name_in, model_part_name_out, forward):
        super().initialize()

        if forward:
            mp_in = model.get_model_part(model_part_name_in)
            n_in = mp_in.size
            n_out = self.n_depth * n_in
            self.n_from = n_in
            self.n_to = n_out

            coords_in = np.column_stack((mp_in.x0, mp_in.y0, mp_in.z0))

            coords_out = np.zeros((n_out, 3))
            ids_out = np.arange(n_out)

            for i_d in range(self.n_depth):  # new nodes ordered per depth coordinate
                i_start = i_d * n_in
                i_end = i_start + n_in

                coords_out[i_start:i_end, self.dir_other[0]] = coords_in[:, self.dir_other[0]]
                coords_out[i_start:i_end, self.dir_other[1]] = coords_in[:, self.dir_other[1]]
                coords_out[i_start:i_end, self.dir_depth] = self.depth_coords[i_d]

            model.create_model_part(model_part_name_out, coords_out[:, 0],
                                    coords_out[:, 1], coords_out[:, 2], ids_out)
        else:
            raise NotImplementedError('Backward Initialization not implemented for MapperDepth2DTo3D.')

    def __call__(self, args_from, args_to):
        super().__call__(args_from, args_to)

        interface_from, mp_name_from, var = args_from
        interface_to, mp_name_to, _ = args_to

        dimensions = variables_dimensions[var]
        data_from = interface_from.get_variable_data(mp_name_from, var)
        if dimensions == 1:
            data_to = np.tile(data_from, (self.n_depth, 1))
        elif dimensions == 3:
            data_to = np.tile(data_from, (self.n_depth, 1))
            data_to[:, self.dir_depth] = 0
        else:
            raise NotImplementedError(
                f'MapperDepth2DTo3D not implemented for variable of dimension {dimensions}')
        interface_to.set_variable_data(mp_name_to, var, data_to)
