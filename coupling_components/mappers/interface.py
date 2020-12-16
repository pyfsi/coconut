from coconut.coupling_components.component import Component
from coconut.coupling_components import tools
from coconut.coupling_components.tools import create_instance


def create(parameters):
    return MapperInterface(parameters)


class MapperInterface(Component):
    def __init__(self, parameters):
        super().__init__()

        self.parameters = parameters
        self.settings = parameters['settings']
        self.mappers = {}

    def initialize(self, interface_from, interface_to):
        super().initialize()

        # loop over ModelParts and create mappers
        for item_from, item_to in zip(interface_from.parameters,
                                      interface_to.parameters):
            mp_from = interface_from.get_model_part(item_from['model_part'])
            mp_to = interface_to.get_model_part(item_to['model_part'])
            mapper = create_instance(self.settings)
            mapper.initialize(mp_from, mp_to)
            mapper.model_part_names = (mp_from.name, mp_to.name)
            self.mappers[mp_from.name + '_to_' + mp_to.name] = mapper

    def finalize(self):
        super().finalize()

        for mapper in self.mappers.values():
            mapper.finalize()

    def __call__(self, interface_from, interface_to):
        # loop over ModelParts and variables and interpolate
        for pair_from, pair_to in zip(interface_from.model_part_variable_pairs,
                                      interface_to.model_part_variable_pairs):
            mapper = self.mappers[pair_from[0] + '_to_' + pair_to[0]]
            mapper((interface_from, *pair_from), (interface_to, *pair_to))

    def output_solution_step(self):
        for mapper in self.mappers.values():
            mapper.output_solution_step()

    def print_components_info(self, pre):
        tools.print_info(pre, "The component ", self.__class__.__name__, " maps the following model parts:")
        pre = tools.update_pre(pre)
        for i, mapper in enumerate(self.mappers.values()):
            name_from, name_to = mapper.model_part_names
            tmp1, tmp2 = ('└─', '  └─') if i == len(self.mappers.values()) - 1 else ('├─', '│ └─')
            tools.print_info(pre, f"{tmp1}ModelPart '{name_from}' to " + f"ModelPart '{name_to}' with the mapper:")
            mapper.print_components_info(pre + tmp2)
