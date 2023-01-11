from coconut.coupling_components.component import Component


class SolverWrapper(Component):
    def __init__(self, parameters):
        super().__init__()

        self.interface_input = None
        self.interface_output = None

        # time
        # noinspection PyUnresolvedReferences
        self.init_time = self.init_time  # created by decorator time_initialize
        self.run_time = 0.0

        # debug
        self.settings = parameters['settings']
        self.debug = self.settings.get('debug', False)  # save copy of input and output files in every iteration

    def get_interface_input(self):
        return self.interface_input.copy()

    def get_interface_output(self):
        return self.interface_output.copy()

    def get_time_allocation(self):
        return {'init_time': self.init_time, 'run_time': self.run_time}
