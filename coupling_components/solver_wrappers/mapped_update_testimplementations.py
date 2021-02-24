from coconut.tools import create_instance
from coconut.coupling_components.component import Component
from coconut import tools
from coconut import data_structure

import numpy as np

def create(parameters):
    return SolverWrapperMapped_update_testimplementations(parameters)

class SolverWrapperMapped_update_testimplementations(Component):
    def __init__(self, parameters):
        super().__init__()

        self.parameters = parameters
        self.settings = parameters["settings"]

        # create mappers
        self.mapper_interface_input = create_instance(self.settings["mapper_interface_input"])
        self.mapper_interface_output = create_instance(self.settings["mapper_interface_output"])

        self.mapper_interface_input_update_to = create_instance(self.settings["mapper_interface_input_update_to"])
        self.mapper_interface_input_update_from = create_instance(self.settings["mapper_interface_input_update_from"])

        self.interface_input_from = None
        self.interface_input_to = None
        self.interface_output_from = None
        self.interface_output_to = None

    def initialize(self, interface_input_from, interface_input_to, interface_output_from, interface_output_to):
        super().initialize()

        # create input mapper
        self.interface_input_from = interface_input_from.copy()
        self.interface_input_to = self.solver_wrapper.get_interface_input()
        self.mapper_interface_input.initialize(self.interface_input_from, self.interface_input_to)

        # create output mapper
        self.interface_output_to = interface_output_to.copy()
        self.interface_output_from = self.solver_wrapper.get_interface_output()
        self.mapper_interface_output.initialize(self.interface_output_from, self.interface_output_to)

        # create mp input_update_to

        for item_input_to in self.interface_input_to.model:
            parameters_to = [{'model_part': 'new_mp_input_to', 'variables': ['displacement']}]
            self.new_interface_input_to = data_structure.Interface(parameters_to, self.interface_input_to.model)

        # create mp_input_udate_from

        for item_input_from in self.interface_input_from.model:
            parameters_from = [{'model_part': 'new_mp_input_from', 'variables': ['displacement']}]
            self.new_interface_input_from = data_structure.Interface(parameters_from, self.interface_input_from.model)

        # create input_update_to mapper
        self.mapper_interface_input_update_to.initialize(self.interface_output_from, self.new_interface_input_to)

        # create input_update_from mapper
        self.mapper_interface_input_update_from.initialize(self.interface_output_to, self.new_interface_input_from)