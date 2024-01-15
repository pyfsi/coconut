from coconut.tools import create_instance
from coconut.coupling_components.component import Component
from coconut import tools


def create(parameters):
    return SolverWrapperMapped(parameters)


class SolverWrapperMapped(Component):
    mapped = True

    @tools.time_initialize
    def __init__(self, parameters):
        super().__init__()

        # read parameters
        self.parameters = parameters
        self.settings = parameters['settings']

        # create solver
        tools.pass_on_parameters(self.settings, self.settings['solver_wrapper']['settings'],
                                 ['timestep_start', 'delta_t', 'save_restart'])
        self.solver_wrapper = create_instance(self.settings['solver_wrapper'])

        # create mappers
        self.mapper_interface_input = create_instance(self.settings['mapper_interface_input'])
        self.mapper_interface_output = create_instance(self.settings['mapper_interface_output'])

        self.interface_input_from = None
        self.interface_input_to = None
        self.interface_output_to = None

        # time
        # noinspection PyUnresolvedReferences
        self.init_time = self.init_time  # created by decorator time_initialize
        self.run_time = 0.0
        self.save_time = 0.0

    @tools.time_initialize
    def initialize(self, interface_input_from, interface_output_to):
        super().initialize()

        self.solver_wrapper.initialize()

        # initialize input mapper
        self.interface_input_from = interface_input_from.copy()
        self.interface_input_to = self.solver_wrapper.get_interface_input()
        self.mapper_interface_input.initialize(self.interface_input_from, self.interface_input_to)

        # initialize output mapper
        self.interface_output_to = interface_output_to.copy()
        interface_output_from = self.solver_wrapper.get_interface_output()
        self.mapper_interface_output.initialize(interface_output_from, self.interface_output_to)

    def initialize_solution_step(self):
        super().initialize_solution_step()

        self.solver_wrapper.initialize_solution_step()

    @tools.time_solve_solution_step
    def solve_solution_step(self, interface_input_from):
        self.interface_input_from = interface_input_from.copy()
        self.mapper_interface_input(self.interface_input_from, self.interface_input_to)
        interface_output_from = self.solver_wrapper.solve_solution_step(self.interface_input_to)
        self.mapper_interface_output(interface_output_from, self.interface_output_to)
        return self.interface_output_to

    def finalize_solution_step(self):
        super().finalize_solution_step()

        self.solver_wrapper.finalize_solution_step()

    @tools.time_save
    def output_solution_step(self):
        super().output_solution_step()

        self.solver_wrapper.output_solution_step()
        self.mapper_interface_input.output_solution_step()
        self.mapper_interface_output.output_solution_step()

    def finalize(self):
        super().finalize()

        self.solver_wrapper.finalize()
        self.mapper_interface_input.finalize()
        self.mapper_interface_output.finalize()

    def get_interface_input(self):
        # does not contain most recent data
        return self.interface_input_from.copy()

    def get_interface_output(self):
        interface_output_from = self.solver_wrapper.get_interface_output()
        self.mapper_interface_output(interface_output_from, self.interface_output_to)
        return self.interface_output_to.copy()

    def get_time_allocation(self):
        time_allocation = {}
        for time_type in ('init_time', 'run_time', 'save_time'):
            total_time = self.__getattribute__(time_type)
            solver_wrapper_time = self.solver_wrapper.get_time_allocation()[time_type]
            mapper_time = total_time - (
                solver_wrapper_time['total'] if isinstance(solver_wrapper_time, dict) else solver_wrapper_time)
            time_allocation[time_type] = {'total': total_time, 'mapper': mapper_time,
                                          'solver_wrapper': solver_wrapper_time}
        return time_allocation

    def print_components_info(self, pre):
        tools.print_info(pre, 'The component ', self.__class__.__name__, ' maps the following solver wrapper:')
        pre = tools.update_pre(pre)
        self.solver_wrapper.print_components_info(pre + '├─')
        tools.print_info(pre, '├─', 'Input mapper:')
        self.mapper_interface_input.print_components_info(pre + '│ └─')
        tools.print_info(pre, '└─', 'Output mapper:')
        self.mapper_interface_output.print_components_info(pre + '  └─')
