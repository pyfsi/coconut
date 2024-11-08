from coconut.coupling_components.component import Component
from coconut import tools


class SolverWrapper(Component):
    # Solver variable check, should be set in subclass
    accepted_out_var = None
    accepted_in_var = None
    error_out = None
    error_in = None

    def __init__(self, parameters):
        super().__init__()

        self.settings = parameters['settings']
        self.type = parameters['type']

        if not self.mapped:
            self.interface_input = None
            self.interface_output = None

            # coupling convergence
            self.check_coupling_convergence = False  # do check of convergence after 1 iteration
            self.coupling_convergence = True  # indicates if solver has converged after 1 iteration
            self.print_coupling_convergence = self.settings.get('print_coupling_convergence', False)

        # time
        # noinspection PyUnresolvedReferences
        self.init_time = self.init_time  # created by decorator time_initialize
        self.run_time = 0.0
        self.save_time = 0.0

        # debug
        self.debug = self.settings.get('debug', False)  # save copy of input and output files in every iteration

    def initialize(self):
        super().initialize()

        # check variables in parameter.json file
        if self.accepted_out_var is not None:
            for mp in self.settings['interface_output']:
                for var in mp['variables']:
                    if var not in self.accepted_out_var:
                        raise NameError(self.error_out)

        if self.accepted_in_var is not None:
            for mp in self.settings['interface_input']:
                for var in mp['variables']:
                    if var not in self.accepted_in_var:
                        raise NameError(self.error_in)

    def initialize_solution_step(self):
        super().initialize_solution_step()

        if not self.mapped:
            self.coupling_convergence = False

    def get_interface_input(self):
        return self.interface_input.copy()

    def get_interface_output(self):
        return self.interface_output.copy()

    def get_time_allocation(self):
        return {'init_time': self.init_time, 'run_time': self.run_time, 'save_time': self.save_time}
