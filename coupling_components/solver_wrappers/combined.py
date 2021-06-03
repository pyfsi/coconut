from coconut.tools import create_instance
from coconut.coupling_components.component import Component
from coconut import tools


def create(parameters):
    return SolverWrapperCombined(parameters)


class SolverWrapperCombined(Component):
    def __init__(self, parameters):
        super().__init__()

        self.parameters = parameters
        self.settings = parameters["settings"]

        # create solvers
        for sol_wrapper_param in self.settings["solver_wrappers"]:
            tools.pass_on_parameters(self.settings, sol_wrapper_param["settings"],
                                     ["timestep_start", "delta_t"])
        nr_mapped_sol_wrappers = 0
        for index, sol_wrapper_param in enumerate(self.settings["solver_wrappers"]):
            if sol_wrapper_param["type"] == "solver_wrappers.mapped":
                nr_mapped_sol_wrappers += 1
            else:
                self.master_sol_index = index

        nr_sol_wrappers = len(self.settings["solver_wrappers"])
        if not nr_mapped_sol_wrappers == nr_sol_wrappers - 1:
            raise RuntimeError(f'Required  number of  mapped solver wrappers: {nr_sol_wrappers - 1}, '
                               f'but {nr_mapped_sol_wrappers} is provided.')

        self.solver_wrapper_list = []
        for sol_wrapper_param in self.settings["solver_wrappers"]:
            self.solver_wrapper_list.append(create_instance(sol_wrapper_param))

        self.mapped_solver_wrapper_list = []
        for index, sol_wrapper in enumerate(self.solver_wrapper_list):
            if not index == self.master_sol_index:
                self.mapped_solver_wrapper_list.append(sol_wrapper)

        self.master_solver_wrapper = self.solver_wrapper_list[self.master_sol_index]

        self.interface_output = None

        # run time
        self.run_time = 0.0

    def initialize(self):
        super().initialize()
        self.master_solver_wrapper.initialize()
        interface_input_from = self.master_solver_wrapper.get_interface_input()
        interface_output_to = self.master_solver_wrapper.get_interface_output().copy()

        for sol_wrapper in self.mapped_solver_wrapper_list:
            sol_wrapper.initialize(interface_input_from, interface_output_to)

    def initialize_solution_step(self):
        super().initialize_solution_step()
        for sol_wrapper in self.solver_wrapper_list:
            sol_wrapper.initialize_solution_step()

    @tools.time_solve_solution_step
    def solve_solution_step(self, interface_input):
        self.interface_output = self.master_solver_wrapper.solve_solution_step(interface_input)
        for sol_wrapper in self.mapped_solver_wrapper_list:
            other_interface_output = sol_wrapper.solve_solution_step(interface_input)
            self.interface_output += other_interface_output

        return self.interface_output

    def finalize_solution_step(self):
        super().finalize_solution_step()
        for sol_wrapper in self.solver_wrapper_list:
            sol_wrapper.finalize_solution_step()

    def finalize(self):
        super().finalize()
        for sol_wrapper in self.solver_wrapper_list:
            sol_wrapper.finalize()

    def output_solution_step(self):
        super().output_solution_step()
        for sol_wrapper in self.solver_wrapper_list:
            sol_wrapper.output_solution_step()

    def get_interface_input(self):
        # does not contain most recent data
        return self.master_solver_wrapper.get_interface_input()

    def get_interface_output(self):
        self.interface_output = self.master_solver_wrapper.get_interface_output()
        for sol_wrapper in self.mapped_solver_wrapper_list:
            other_interface_output = sol_wrapper.get_interface_output()
            self.interface_output += other_interface_output
        return self.interface_output

    def print_components_info(self, pre):
        tools.print_info(pre, "The component ", self.__class__.__name__, " combines the following solver wrappers:")
        pre = tools.update_pre(pre)
        tools.print_info(pre, '├─', "Mapped solver wrappers:")
        for sol_wrapper in self.mapped_solver_wrapper_list[:-1]:
            sol_wrapper.print_components_info(pre + '│ ├─')
        self.solver_wrapper_list[-1].print_components_info(pre + '│ └─')

        tools.print_info(pre, '└─', "Master solver wrapper:")
        self.master_solver_wrapper.print_components_info(pre + '  └─')
