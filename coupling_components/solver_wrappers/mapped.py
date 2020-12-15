from coconut.coupling_components.tools import create_instance
from coconut.coupling_components.component import Component
from coconut.coupling_components import tools

""" proposed changes to mapped.py
- do initialization of mappers in Initialize method, would be more logical
- remove all set_interface_input/Output methods?
- use copy in get_interface_input/Output methods?
    and just refer to actual solver wrapper in SolverWrapperMapped
- all Interfaces are stored in this mapper, e.g. self.interface_output_to and 3 others;
    I see no reason for this; furthermore, it is only useful to store it if you take copies all the time
- output_solution_step is barely used; what's the deal with it??
"""


def create(parameters):
    return SolverWrapperMapped(parameters)


class SolverWrapperMapped(Component):
    def __init__(self, parameters):
        super().__init__()

        # Read parameters
        self.parameters = parameters
        self.settings = parameters["settings"]

        # Create solver
        self.solver_wrapper = create_instance(self.settings["solver_wrapper"])

        # run time
        self.run_time = 0.0

    def initialize(self):
        super().initialize()

        self.solver_wrapper.initialize()

    def initialize_solution_step(self):
        super().initialize_solution_step()

        self.solver_wrapper.initialize_solution_step()

    @tools.time_solve_solution_step
    def solve_solution_step(self, interface_input_from):
        self.interface_input_from = interface_input_from
        self.mapper_interface_input(self.interface_input_from, self.interface_input_to)
        self.interface_output_from = self.solver_wrapper.solve_solution_step(self.interface_input_to)
        self.mapper_interface_output(self.interface_output_from, self.interface_output_to)
        return self.interface_output_to

    def finalize_solution_step(self):
        super().finalize_solution_step()

        self.solver_wrapper.finalize_solution_step()

    def finalize(self):
        super().finalize()

        self.solver_wrapper.finalize()
        self.mapper_interface_input.finalize()
        self.mapper_interface_output.finalize()

    def output_solution_step(self):
        super().output_solution_step()

        self.solver_wrapper.output_solution_step()
        self.mapper_interface_input.output_solution_step()
        self.mapper_interface_output.output_solution_step()

    def get_interface_input(self):
        # Does not contain most recent data
        # *** shouldn't this just call the underlying solver wrapper?
        return self.interface_input_from

    def set_interface_input(self, interface_input_from):
        # Create input mapper
        self.interface_input_from = interface_input_from.deepcopy()
        self.interface_input_to = self.solver_wrapper.get_interface_input()

        self.mapper_interface_input = create_instance(self.settings["mapper_interface_input"])
        self.mapper_interface_input.initialize(self.interface_input_from, self.interface_input_to)

    def get_interface_output(self):
        self.interface_output_from = self.solver_wrapper.get_interface_output()
        self.mapper_interface_output(self.interface_output_from, self.interface_output_to)
        return self.interface_output_to.deepcopy()

    def set_interface_output(self, interface_output_to):
        # Create output mapper
        self.interface_output_to = interface_output_to.deepcopy()
        self.interface_output_from = self.solver_wrapper.get_interface_output()

        self.mapper_interface_output = create_instance(self.settings["mapper_interface_output"])
        self.mapper_interface_output.initialize(self.interface_output_from, self.interface_output_to)

    def print_components_info(self, pre):
        tools.print_info(pre, "The component ", self.__class__.__name__, " maps the following solver wrapper:")
        pre = tools.update_pre(pre)
        self.solver_wrapper.print_components_info(pre + '├─')
        tools.print_info(pre, '├─', "Input mapper:")
        self.mapper_interface_input.print_components_info(pre + '│ └─')
        tools.print_info(pre, '└─', "Output mapper:")
        self.mapper_interface_output.print_components_info(pre + '  └─')
