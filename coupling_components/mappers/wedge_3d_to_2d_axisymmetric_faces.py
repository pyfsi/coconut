from coconut.coupling_components.mappers.transformer import MapperTransformer
from coconut.data_structure import variables_dimensions
from scipy.spatial import cKDTree

import numpy as np


def create(parameters):
    return MapperWedge3DToAxisymmetric2D(parameters)

# TODO: mention in docs that these mappers cannot handle singular points (i.e. with r = 0, e.g. balloon)

class MapperWedge3DToAxisymmetric2D(MapperTransformer):
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
        super().initialize()

        if forward:

            mp_in = model.get_model_part(model_part_name_in)
            n_in = mp_in.size
            n_out = n_in
            self.n_from = n_in
            self.n_to = n_out

            self.coords_in = np.column_stack((mp_in.x0, mp_in.y0, mp_in.z0))

            coords_tmp = np.zeros((n_out, 3))
            coords_out = np.zeros((n_out, 3))
            ids_out = np.arange(n_out)

            self.nearest = np.zeros((self.n_to,)).astype(int)

            j = 0
            i_to = 0
            for i in range(n_in):
                coords_tmp[j,: ] = self.coords_in[i,:]
                self.nearest[i] = i_to
                j += 1
                i_to += 1

            coords_out[:, 0] = coords_tmp[:, 0]
            coords_out[:, 1] = coords_tmp[:, 1]/np.cos(self.angle)
            coords_out[:, 2] = 0

            model.create_model_part(model_part_name_out, coords_out[:, 0],
                                    coords_out[:, 1], coords_out[:, 2], ids_out)
        else:
            raise NotImplementedError('Forward Initialization not implemented for MapperWedge3DTo2DAxisymmetric.')

    def __call__(self, args_from, args_to):
        super().__call__(args_from, args_to)

        interface_from, mp_name_from, var = args_from
        interface_to, mp_name_to, _ = args_to

        # print("args_from_faces")
        # # print(interface_from.parameters)
        # for i in interface_from.parameters:
        #     mp_from = interface_from.get_model_part(i["model_part"])
        #     print(mp_from.y0)

        dimensions = variables_dimensions[var]
        data_from = interface_from.get_variable_data(mp_name_from, var)
        if dimensions == 1:
            data_to = np.tile(data_from,(1,1))
            for i_to in range(len(data_to)):
                data_to[i_to] = data_from[self.nearest[i_to]]

        elif dimensions == 3:
            data_to = np.tile(data_from, (1,1))
            for i_to in range(len(data_to)):
                data_to[i_to] = data_from[self.nearest[i_to]]

            r = data_to[:, self.dir_r].copy()
            data_to[:, self.dir_r] = r / np.cos(self.angle)
            data_to[:, self.dir_3d] = 0


        else:
            raise NotImplementedError(f'MapperWedge3DTo2DAxisymmetric not implemented for variable of dimension {dimensions}')
        interface_to.set_variable_data(mp_name_to, var, data_to)
