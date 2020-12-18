from coconut.coupling_components.component import Component
from coconut import data_structure


def create(parameters):
    raise NotImplementedError('this class can only be used as super-class')


class MapperTransformer(Component):
    def __init__(self, parameters):
        super().__init__()

        self.settings = parameters['settings']
        self.interpolator = False

    def __call__(self, args_from, args_to):
        # check variables
        var_from, var_to = args_from[2], args_to[2]
        if var_from not in data_structure.variables_dimensions:
            raise NameError(f'variable "{var_from}" does not exist')
        if var_from != var_to:
            raise TypeError('variables must be equal')
