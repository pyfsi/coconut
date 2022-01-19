from coconut.coupling_components.mappers.transformer import MapperTransformer

import numpy as np

def create(parameters):
    return Mapper_offset_displacement_velocity(parameters)


class Mapper_offset_displacement_velocity(MapperTransformer):
    def __init__(self, parameters):
        super().__init__(parameters)


        self.depth = self.settings['initial_depth']
        if type(self.depth) != float:
            raise TypeError('coords_min must be a float')

    def initialize(self, model, model_part_name_in, model_part_name_out, forward):
        if forward:
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
            raise NotImplementedError('Forward Initialization not implemented for initial_load')

    def __call__(self, args_from, args_to):
        super().__call__(args_from, args_to)
        interface_from, mp_name_from, var = args_from
        interface_to, mp_name_to, _ = args_to
        data_from = interface_from.get_variable_data(mp_name_from, var)

        for input_from in interface_from.parameters:
            mp_from = interface_from.get_model_part(input_from['model_part'])

        data_to = np.zeros((3,len(mp_from.y0)))
        data_to[0] = data_from[:,0] * 0.0000000000001
        data_to[1] = data_from[:,1]
        for i in range((len(mp_from.z0))):
            if mp_from.z0[i] < 0:
                data_to[2,i] = mp_from.z0[i] + self.depth
            else:
                data_to[2,i] = mp_from.z0[i] - self.depth
        # for i in range((len(mp_from.z0))):
        #     if mp_from.z0[i] < 0:
        #         data_to[2,i] = -data_from[2,i]
        #     else:
        #         data_to[2,i] = data_from[2,i]
        data_to = np.transpose(data_to)
        # print("data_aftermapper")
        # print(data_to)
        interface_to.set_variable_data(mp_name_to, var, data_to)


