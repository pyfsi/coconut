from coconut.coupling_components import tools
from coconut.coupling_components.component import Component
from coconut.coupling_components.tools import create_instance
from coconut import data_structure
from coconut.data_structure import variables_dimensions


def create(parameters):
    return MapperCombined(parameters)


class MapperCombined(Component):
    def __init__(self, parameters):
        """
        This mapper combines 1 interpolator mapper
        with 0 or more transformer mappers, in
        arbitrary order.

        If more than 1 mapper is present, inner
        ModelParts are created and stored in the
        attribute model_parts.
        """
        super().__init__()

        self.settings = parameters['settings']

        # create all mappers
        self.mappers = []
        for par in self.settings['mappers']:
            self.mappers.append(create_instance(par))

        # initialization
        self.interfaces = None

        # check that exactly one mapper is an interpolator
        counter = 0
        for i, mapper in enumerate(self.mappers):
            if mapper.interpolator:
                self.index = i
                counter += 1
        if counter != 1:
            raise ValueError(f'{counter} interpolators found instead of 1')

    def initialize(self, model_part_from, model_part_to):
        super().initialize()

        self.model = data_structure.Model()
        n = len(self.mappers)

        # initialize upstream transformers
        mp = model_part_from
        mp_names_from = [mp.name]
        self.model.create_model_part(mp.name, mp.x0, mp.y0, mp.z0, mp.id)
        for i in range(self.index):
            mp_names_from.append(f'mp_from_{i + 1}')
            self.mappers[i].initialize(self.model, mp_names_from[i],
                                       mp_names_from[i + 1], forward=True)

        # initialize downstream transformers
        mp = model_part_to
        mp_names_to = [mp.name]
        self.model.create_model_part(mp.name, mp.x0, mp.y0, mp.z0, mp.id)
        for i in range(n - 1 - self.index):
            mp_names_to.append(f'mp_to_{i + 1}')
            self.mappers[n - 1 - i].initialize(self.model, mp_names_to[i],
                                               mp_names_to[i + 1], forward=False)

        # initialize interpolator
        mp_from = self.model.get_model_part(mp_names_from[-1])
        mp_to = self.model.get_model_part(mp_names_to[-1])
        self.mappers[self.index].initialize(mp_from, mp_to)

        self.mp_names = mp_names_from + mp_names_to[::-1]

    def __call__(self, args_from, args_to):
        # unpack arguments
        interface_from, mp_name_from, var_from = args_from
        interface_to, mp_name_to, var_to = args_to

        # check variables
        if var_from not in variables_dimensions:
            raise NameError(f'Variable "{var_from}" does not exist')
        if var_from != var_to:
            raise TypeError('Variables must be equal')
        var = var_from

        # create Interfaces
        interfaces = []
        for mp_name in self.mp_names:
            parameters = [{'model_part': mp_name, 'variables': [var]}]
            interfaces.append(data_structure.Interface(parameters, self.model))

        # map data
        data_from = interface_from.get_variable_data(mp_name_from, var)
        interfaces[0].set_variable_data(self.mp_names[0], var, data_from)

        for i, mapper in enumerate(self.mappers):
            mapper((interfaces[i], self.mp_names[i], var),
                   (interfaces[i + 1], self.mp_names[i + 1], var))

        data_to = interfaces[-1].get_variable_data(self.mp_names[-1], var)
        interface_to.set_variable_data(mp_name_to, var, data_to)

    def output_solution_step(self):
        for mapper in self.mappers:
            mapper.output_solution_step()

    def print_components_info(self, pre):
        tools.print_info(pre, "The component ", self.__class__.__name__, " combining the following mappers:")
        tools.print_components_info(pre, self.mappers)
