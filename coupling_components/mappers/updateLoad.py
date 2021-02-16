from coconut.coupling_components.mappers.transformer import MapperTransformer
from coconut.data_structure import variables_dimensions
from scipy.spatial import cKDTree
from coconut import data_structure
import numpy as np
from scipy import interpolate

import numpy as np

def create(parameters):
    return MapperLoadUpDate(parameters)


class MapperLoadUpDate(MapperTransformer):
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



    def initialize(self, model, model_part_name_in, model_part_name_out, forward):
        if not forward:
            super().initialize()

            self.mp_input_to = model.get_model_part(model_part_name_in)

            #TODO: Here should come a call from the interface_input_from, but to test it, is it at the moment not possible.
            model0 = data_structure.Model()
            x0 = np.linspace(0, 2, 5)
            y0 = np.ones(5)
            z0 = np.zeros(5)
            ids0 = np.arange(5)
            self.mp_input_from = model0.create_model_part('mp_input_from', x0, y0, z0, ids0)
            self.v_min = min(self.mp_input_from.x0)
            self.v_max = max(self.mp_input_from.x0)
            parameters = [{'model_part': 'mp_input_from', 'variables': ['pressure']}]
            self.interface_input_from = data_structure.Interface(parameters,model0)

            for input_from in self.interface_input_from.parameters:
                tmp = self.interface_input_from.get_model_part(input_from['model_part'])
            v_min = min(tmp.x0)
            v_max = max(tmp.x0)

            n_in = self.mp_input_to.size
            n_out = n_in
            self.n_from = n_out
            self.n_to = n_in

            self.coords_in = np.column_stack((self.mp_input_to.x0, self.mp_input_to.y0, self.mp_input_to.z0))
            self.x_orig = self.mp_input_to.x0.flatten()
            self.y_orig = self.mp_input_to.y0.flatten()
            self.f = interpolate.interp1d(self.x_orig,self.y_orig)

            # self.tree= cKDTree(self.coords_in, balanced_tree = False)
            # coords = self.coords_in.copy()

            coords_out = np.zeros((n_out, 3))
            ids_out = np.arange(n_out)

            # self.nearest = np.zeros((self.n_to,)).astype(int)

            coords_out[:, 0] = np.linspace(v_min, v_max,self.mp_input_to.size)
            coords_out[:, 1] =self.f(coords_out[:,0])
            coords_out[:, 2] = 0

            self.mp_out = model.create_model_part(model_part_name_out, coords_out[:, 0],
                                    coords_out[:, 1], coords_out[:, 2], ids_out)

        else:
            raise NotImplementedError('Backward Initialization not implemented for LoadUpdate')

    def __call__(self, args_from, args_to):
        super().__call__(args_from, args_to)

        interface_from, mp_name_from, var = args_from
        interface_to, mp_name_to, _ = args_to

        dimensions = variables_dimensions[var]
        self.data_from = interface_from.get_variable_data(mp_name_from, var)

        if dimensions == 1:
            data = self.data_from.flatten()
            g = interpolate.interp1d(self.mp_out.x0, data)
            data_to = np.zeros((self.mp_input_to.size, 1))
            for i in range(self.mp_input_to.size):
                if self.mp_input_to.x0[i] < self.v_min or self.mp_input_to.x0[i] >self.v_max:
                    data_to[i] = 0
                else:
                    data_to[i] = g(self.mp_input_to.x0[i])

        elif dimensions == 3:
            data_tree = self.data_from
            data_x = data_tree[ :,0].flatten()
            data_y = data_tree[ :,1].flatten()
            data_z = data_tree[ :,2].flatten()
            e1 = interpolate.interp1d(self.mp_out.x0, data_x)
            e2 = interpolate.interp1d(self.mp_out.x0, data_y)
            e3 = interpolate.interp1d(self.mp_out.x0, data_z)
            data_to = np.zeros((self.mp_input_to.size, 3))
            for j in range(self. mp_input_to.size):
                if self.mp_input_to.x0[j] < self.v_min or self.mp_input_to.x0[j] >self.v_max:
                    data_to[j,0] = 0
                    data_to[j,1] = 0
                    data_to[j,2] = 0
                else:
                    data_to[j,0] = e1(self.mp_input_to.x0[j])
                    data_to[j, 1] = e2(self.mp_input_to.x0[j])
                    data_to[j, 2] = e3(self.mp_input_to.x0[j])


        else:
            raise NotImplementedError(
                f'MapperUpdateLoad not implemented for variable of dimension {dimensions}')
        interface_to.set_variable_data(mp_name_to, var, data_to)