from coconut.coupling_components.component import Component
from coconut import tools


class SolverWrapper(Component):

    def __init__(self, parameters):
        super().__init__()

        self.settings = parameters['settings']
        self.type = parameters["type"]

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
        warning = None
        skip = False
        #TODO: merge variable control from pc & cht solver wrappers to here
        if "mapped" in self.type or "pc" in self.type or "cht" in self.type:
            skip = True
        elif "fluent" in self.type or "openfoam" in self.type:
            accepted_out_var = ["pressure", "traction"]
            accepted_in_var = ["displacement"]
            error_out = f"Only permitted output variables are pressure and traction for {self.type}."
            error_in = f"Only permitted input variable is displacement for {self.type}."
        elif "abaqus" in self.type or "kratos_structure" in self.type:
            accepted_out_var = ["displacement"]
            accepted_in_var = ["pressure", "traction"]
            error_out = f"Only permitted output variable is displacement for {self.type}."
            error_in = f"Only permitted input variables are pressure and traction for {self.type}."
        else:
            warning = "Variables could not be checked as solver_wrapper was not recognized."

        if warning is None and skip is False:
            for mp in self.settings['interface_output']:
                for var in mp['variables']:
                    if var not in accepted_out_var:
                        raise NameError(error_out)

            for mp in self.settings['interface_input']:
                for var in mp['variables']:
                    if var not in accepted_in_var:
                        raise NameError(error_in)
        elif warning is not None and skip is False:
            tools.print_info(warning, layout='warning')

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
