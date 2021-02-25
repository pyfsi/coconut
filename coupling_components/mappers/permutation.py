from coconut.coupling_components.mappers.transformer import MapperTransformer
from coconut.data_structure import variables_dimensions

import numpy as np


def create(parameters):
    return MapperPermutation(parameters)


class MapperPermutation(MapperTransformer):
    def __init__(self, parameters):
        """ TODO: move all this info to docs
        This is not an interpolator, but a transformer.
        This is denoted by setting the self.interpolator
        attribute to False in the __init__.

        The difference is that a transformer is
        initialized from only one side (from or to)
        and returns the ModelPart corresponding
        to the forward or backward transformation.
        An interpolator is initialized from both
        sides and returns nothing.

        It can be initialized from both sides,
        based on the forward parameter.
        If forward == True, then the model_part_from
        is expected as input.
        Else, the model_part_to is expected.

        The historical variables are not changed,
        simply copied from input to output ModelPart.
        The coordinates are permutated according to the
        permutation parameter (list of ints) in the
        JSON file.
        """
        super().__init__(parameters)

        self.permutation = self.settings['permutation']
        if type(self.permutation) != list:
            raise TypeError('Parameter permutation must be a list')
        if set(self.permutation) != {0, 1, 2}:
            raise ValueError('Parameter permutation must be a permutation of [0, 1, 2]')

    def initialize(self, model, model_part_name_in, model_part_name_out, forward):
        super().initialize()

        permutation = self.permutation if forward else np.argsort(self.permutation)

        mp_in = model.get_model_part(model_part_name_in)
        coords_in = np.column_stack((mp_in.x0, mp_in.y0, mp_in.z0))
        coords_out = coords_in[:, permutation]
        model.create_model_part(model_part_name_out, coords_out[:, 0],
                                coords_out[:, 1], coords_out[:, 2], np.arange(mp_in.size))

    def __call__(self, args_from, args_to):
        super().__call__(args_from, args_to)

        interface_from, mp_name_from, var = args_from
        interface_to, mp_name_to, _ = args_to

        dimensions = variables_dimensions[var]
        data = interface_from.get_variable_data(mp_name_from, var)
        if dimensions == 1:
            interface_to.set_variable_data(mp_name_to, var, data)
        elif dimensions == 3:
            interface_to.set_variable_data(mp_name_to, var, data[:, self.permutation])
        else:
            raise NotImplementedError(f'Permutation not implemented for variable of dimension {dimensions}')
