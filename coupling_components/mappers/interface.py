from coconut.coupling_components.component import Component
from coconut.coupling_components import tools
from coconut.coupling_components.tools import create_instance


def create(parameters):
    return MapperInterface(parameters)


# Class MapperInterface: Interface interpolation with same interpolator type for all ModelParts and variables.
class MapperInterface(Component):
    def __init__(self, parameters):
        super().__init__()

        self.parameters = parameters
        self.settings = parameters['settings']
        self.mappers = {}

    def Initialize(self, interface_from, interface_to):
        super().Initialize()

        # loop over ModelParts and create mappers
        for item_from, item_to in zip(interface_from.parameters,
                                      interface_to.parameters):
            mp_from = interface_from.get_model_part(item_from['model_part'])
            mp_to = interface_to.get_model_part(item_to['model_part'])
            mapper = create_instance(self.settings)
            mapper.Initialize(mp_from, mp_to)
            mapper.model_part_names = (mp_from.name, mp_to.name)
            self.mappers[mp_from.name + '_to_' + mp_to.name] = mapper

    def Finalize(self):
        super().Finalize()

        for mapper in self.mappers.values():
            mapper.Finalize()

    def __call__(self, interface_from, interface_to):
        # loop over ModelParts and variables and interpolate
        for pair_from, pair_to in zip(interface_from.model_part_variable_pairs,
                                      interface_to.model_part_variable_pairs):
            mapper = self.mappers[pair_from[0] + '_to_' + pair_to[0]]
            mapper((interface_from, *pair_from), (interface_to, *pair_to))

    def OutputSolutionStep(self):
        for mapper in self.mappers:
            mapper.OutputSolutionStep()

    def PrintInfo(self, pre):
        tools.Print(pre, "The component ", self.__class__.__name__, " maps the following model parts:")
        pre = tools.UpdatePre(pre)
        for i, mapper in enumerate(self.mappers.values()):
            name_from, name_to = mapper.model_part_names
            tmp1, tmp2 = ('└─', '  └─') if i == len(self.mappers.values()) - 1 else ('├─', '│ └─')
            tools.Print(pre, f"{tmp1}ModelPart '{name_from}' to " + f"ModelPart '{name_to}' with the mapper:")
            mapper.PrintInfo(pre + tmp2)


        #     tools.Print(pre, f"├─ModelPart '{mapper.mp_from.name}' to " +
        #                 f"ModelPart '{mapper.mp_to.name}' with the mapper:")
        #     mapper.PrintInfo(pre + '│ └─')
        # tools.Print(pre, f"└─ModelPart '{self.mappers[-1].mp_from.name}' to " +
        #             f"ModelPart '{self.mappers[-1].mp_to.name}' with the mapper:")
        # self.mappers[-1].PrintInfo(pre + '  └─')
