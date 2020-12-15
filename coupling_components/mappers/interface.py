from coconut import data_structure
from coconut.coupling_components.component import Component
from coconut.coupling_components import tools
from coconut.coupling_components.tools import CreateInstance


def Create(parameters):
    return MapperInterface(parameters)


# Class MapperInterface: Interface interpolation with same interpolator type for all modelparts and variables.
class MapperInterface(Component):
    def __init__(self, parameters):
        super().__init__()

        self.parameters = parameters
        self.settings = parameters["settings"]
        self.mappers = []

    def initialize(self, interface_from, interface_to):
        super().initialize()

        # Loop over ModelParts and create mappers
        self.mappers = []
        self.keys = []
        for item_from, item_to in zip(interface_from.model_parts_variables,
                                      interface_to.model_parts_variables):
            key_from = item_from[0]
            key_to = item_to[0]
            self.keys.append((key_from, key_to))  # for PrintInfo

            self.mappers.append(CreateInstance(self.settings))
            self.mappers[-1].initialize(interface_from.model[key_from],
                                        interface_to.model[key_to])

    def finalize(self):
        super().finalize()

        for mapper in self.mappers:
            mapper.finalize()

    def __call__(self, interface_from, interface_to):
        # Loop over ModelParts and Variables to interpolate
        for i, mapper in enumerate(self.mappers):
            key_from, variables_from = interface_from.model_parts_variables[i]
            key_to, variables_to = interface_to.model_parts_variables[i]
            model_part_from = interface_from.model[key_from]
            model_part_to = interface_to.model[key_to]
            for var_from, var_to in zip(variables_from.list(), variables_to.list()):
                mapper((model_part_from, vars(data_structure)[var_from.GetString()]),
                       (model_part_to, vars(data_structure)[var_to.GetString()]))

    def output_solution_step(self):
        for mapper in self.mappers:
            mapper.output_solution_step()

    def print_components_info(self, pre):
        tools.print_info(pre, "The component ", self.__class__.__name__, " maps the following model parts:")
        pre = tools.update_pre(pre)
        for i, mapper in enumerate(self.mappers[:-1]):
            tools.print_info(pre, f"├─ModelPart '{self.keys[i][0]}' to ModelPart '{self.keys[i][1]}' with the mapper:")
            mapper.print_components_info(pre + '│ └─')
        tools.print_info(pre, f"└─ModelPart '{self.keys[-1][0]}' to ModelPart '{self.keys[-1][1]}' with the mapper:")
        self.mappers[-1].print_components_info(pre + '  └─')
