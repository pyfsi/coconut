from coconut.coupling_components.mappers.transformer import MapperTransformer
from coconut.data_structure import variables_dimensions
from scipy.spatial import cKDTree

import numpy as np


def create(parameters):
    return MapperAxisymmetric2DToWedge3D(parameters)

# TODO: mention in docs that these mappers cannot handle singular points (i.e. with r = 0, e.g. balloon)

class MapperAxisymmetric2DToWedge3D(MapperTransformer):
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
        self.angle = np.radians(2.5)


    def initialize(self, model, model_part_name_in, model_part_name_out, forward):
        if not forward:
            super().initialize()

            mp_in = model.get_model_part(model_part_name_in)
            n_in = mp_in.size
            n_out = n_in // 2
            self.n_from = n_out
            self.n_to = n_in

            self.coords_in = np.column_stack((mp_in.x0, mp_in.y0, mp_in.z0))
            self.tree= cKDTree(self.coords_in, balanced_tree = False)
            coords = self.coords_in.copy()


            coords_tmp = np.zeros((n_out, 3))
            coords_out = np.zeros((n_out, 3))
            ids_out = np.arange(n_out)

            self.nearest = np.zeros((self.n_to,)).astype(int)

            j = 0
            i_from = 0
            for i in range(n_in):
                if coords[i, 2] > 0:
                    coords_tmp[j, :] = self.coords_in[i, :]

                    j += 1
                    self.nearest[i] = i_from
                    coords[i, 2] = -coords[i, 2]
                    dd, ii = self.tree.query(coords[i], k=1)
                    self.nearest[ii] = i_from
                    i_from += 1

            coords_out[:, 0] = coords_tmp[:, 0]
            coords_out[:, 1] = np.cos(self.angle) * coords_tmp[:, 1] + np.sin(self.angle) * coords_tmp[:, 2]
            coords_out[:, 2] = 0

            model.create_model_part(model_part_name_out, coords_out[:, 0],
                                    coords_out[:, 1], coords_out[:, 2], ids_out)
        else:
            raise NotImplementedError('Backward Initialization not implemented for Mapper2DAxisymmetricToWedge3D.')

    def __call__(self, args_from, args_to):
        super().__call__(args_from, args_to)

        interface_from, mp_name_from, var = args_from
        interface_to, mp_name_to, _ = args_to

        dimensions = variables_dimensions[var]
        data_from = interface_from.get_variable_data(mp_name_from, var)

        if dimensions == 1:
            data_to = np.tile(data_from,(2,1))

            for i_to in range(len(data_to)):
                data_to[i_to] = data_from[self.nearest[i_to]]

        elif dimensions == 3:
            data_to = np.tile(data_from, (2, 1))

            for i_to in range(len(data_to)):
                data_to[i_to] = data_from[self.nearest[i_to]]

            z0 = self.coords_in[:, self.dir_3d]
            r = data_to[:, self.dir_r].copy()
            z = data_to[:, self.dir_3d].copy()

            data_to[:, self.dir_r] = r * np.cos(self.angle)

            for i in range(self.n_to):
                if z0[i]<0:
                    z[i] = -r[i]*np.sin(self.angle)
                else:
                    z[i] = r[i] * np.sin(self.angle)

            data_to[:, self.dir_3d] = z

        else:
            raise NotImplementedError(f'Mapper2DAxisymmetricToWedge3D not implemented for variable of dimension {dimensions}')
        interface_to.set_variable_data(mp_name_to, var, data_to)
