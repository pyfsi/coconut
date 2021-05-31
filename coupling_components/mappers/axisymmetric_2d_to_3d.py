from coconut.coupling_components.mappers.transformer import MapperTransformer
from coconut.data_structure import variables_dimensions

import numpy as np


def create(parameters):
    return MapperAxisymmetric2DTo3D(parameters)

# TODO: mention in docs that these mappers cannot handle singular points (i.e. with r = 0, e.g. balloon)

class MapperAxisymmetric2DTo3D(MapperTransformer):
    def __init__(self, parameters):
        super().__init__(parameters)

        # get axial and radial directions
        dirs = ['x', 'y', 'z']
        if self.settings['direction_axial'] not in dirs:
            raise ValueError(f'invalid axial_direction {self.settings["direction_axial"]}')
        if self.settings['direction_axial'] not in dirs:
            raise ValueError(f'invalid radial_direction {self.settings["direction_radial"]}')
        self.dir_a = dirs.index(self.settings['direction_axial'])
        self.dir_r = dirs.index(self.settings['direction_radial'])
        self.dir_3d = ({0, 1, 2} - {self.dir_a, self.dir_r}).pop()
        self.angle  = self.settings.get('angle',360) #angle is set in degrees

        # get number of nodes in tangential direction
        self.n_t = self.settings['n_tangential']
        limit = self.angle // 72 + 2
        if self.n_t < limit:
            raise ValueError('minimum value for n_tangential is ' + str(limit))

    def initialize(self, model, model_part_name_in, model_part_name_out, forward):
        super().initialize()

        if forward:
            mp_in = model.get_model_part(model_part_name_in)
            n_in = mp_in.size
            n_out = self.n_t * n_in
            self.n_from = n_in
            self.n_to = n_out

            coords_in = np.column_stack((mp_in.x0, mp_in.y0, mp_in.z0))

            coords_out = np.zeros((n_out, 3))
            ids_out = np.arange(n_out)
            self.theta = np.zeros(n_out)

            for i_t in range(self.n_t):  # new nodes ordered per theta
                if self.angle == 360:
                    # theta =  -np.radians(self.angle / 2) + i_t*np.radians(self.angle)/(self.n_t)
                    theta = i_t*np.radians(self.angle)/(self.n_t)
                else:
                    theta = -np.radians(self.angle / 2) + i_t*np.radians(self.angle)/(self.n_t - 1)
                i_start = i_t * n_in
                i_end = (i_t + 1) * n_in

                coords_out[i_start: i_end, self.dir_a] = coords_in[:, self.dir_a]
                coords_out[i_start: i_end, self.dir_r] = np.cos(theta) * coords_in[:, self.dir_r]
                coords_out[i_start: i_end, self.dir_3d] = np.sin(theta) * coords_in[:, self.dir_r]

                self.theta[i_start: i_end] = theta

            model.create_model_part(model_part_name_out, coords_out[:, 0],
                                    coords_out[:, 1], coords_out[:, 2], ids_out)
        else:
            raise NotImplementedError('Backward Initialization not implemented for MapperAxisymmetric2DTo3D.')

    def __call__(self, args_from, args_to):
        super().__call__(args_from, args_to)

        interface_from, mp_name_from, var = args_from
        interface_to, mp_name_to, _ = args_to

        dimensions = variables_dimensions[var]
        data_from = interface_from.get_variable_data(mp_name_from, var)
        if dimensions == 1:
            data_to = np.tile(data_from, (self.n_t, 1))
        elif dimensions == 3:
            data_to = np.tile(data_from, (self.n_t, 1))
            r = data_to[:, self.dir_r].copy()
            data_to[:, self.dir_r] = r * np.cos(self.theta)
            data_to[:, self.dir_3d] = r * np.sin(self.theta)
        else:
            raise NotImplementedError(f'MapperAxisymmetric2DTo3D not implemented for variable of dimension {dimensions}')
        interface_to.set_variable_data(mp_name_to, var, data_to)
