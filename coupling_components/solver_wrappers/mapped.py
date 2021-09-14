from coconut.tools import create_instance
from coconut.coupling_components.component import Component
from coconut import tools


def create(parameters):
    return SolverWrapperMapped(parameters)


class SolverWrapperMapped(Component):
    @tools.time_initialize
    def __init__(self, parameters):
        super().__init__()

        self.parameters = parameters
        self.settings = parameters['settings']

        # create solver
        tools.pass_on_parameters(self.settings, self.settings['solver_wrapper']['settings'],
                                 ['timestep_start', 'delta_t', 'save_restart'])
        self.solver_wrapper = create_instance(self.settings['solver_wrapper'])

        # create mappers
        self.mapper_interface_input = create_instance(self.settings['mapper_interface_input'])
        self.mapper_interface_output = create_instance(self.settings['mapper_interface_output'])
        # self.mapper_interface_output_it = create_instance(self.settings['mapper_interface_output'])

        self.interface_input_from = None
        self.interface_input_to = None
        self.interface_output_to = None

        # time
        self.init_time = 0.0
        self.run_time = 0.0
        self.iteration = 0

    @tools.time_initialize
    def initialize(self, interface_input_from, interface_output_to):
        super().initialize()

        self.solver_wrapper.initialize()

        # create input mapper
        self.interface_input_from = interface_input_from.copy()
        self.interface_input_to = self.solver_wrapper.get_interface_input()
        self.mapper_interface_input.initialize(self.interface_input_from, self.interface_input_to)

        # create output mapper
        self.interface_output_to = interface_output_to.copy()
        interface_output_from = self.solver_wrapper.get_interface_output()
        self.mapper_interface_output.initialize(interface_output_from, self.interface_output_to)
        print("interface_output_initialize")
        print(interface_output_from)

    def initialize_solution_step(self):
        super().initialize_solution_step()

        self.solver_wrapper.initialize_solution_step()

    @tools.time_solve_solution_step
    def solve_solution_step(self, interface_input_from):
        self.iteration +=1
        self.interface_input_from = interface_input_from.copy()
        self.mapper_interface_input(self.interface_input_from, self.interface_input_to)
        interface_output_from = self.solver_wrapper.solve_solution_step(self.interface_input_to)
        self.mapper_interface_output_it = create_instance(self.settings['mapper_interface_output'])
        self.mapper_interface_output_it.initialize(interface_output_from, self.interface_output_to)
        print("interface_output_mapped")
        print(interface_output_from)
        self.mapper_interface_output_it(interface_output_from, self.interface_output_to)
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
        # does not contain most recent data
        return self.interface_input_from

    def get_interface_output(self):
        interface_output_from = self.solver_wrapper.get_interface_output()
        self.mapper_interface_output(interface_output_from, self.interface_output_to)
        return self.interface_output_to

    def print_components_info(self, pre):
        tools.print_info(pre, 'The component ', self.__class__.__name__, ' maps the following solver wrapper:')
        pre = tools.update_pre(pre)
        self.solver_wrapper.print_components_info(pre + '├─')
        tools.print_info(pre, '├─', 'Input mapper:')
        self.mapper_interface_input.print_components_info(pre + '│ └─')
        tools.print_info(pre, '└─', 'Output mapper:')
        self.mapper_interface_output.print_components_info(pre + '  └─')
