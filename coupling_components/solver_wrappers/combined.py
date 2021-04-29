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
        master_sol_index = self.settings["master_solver_index"]
        self.master_solver_wrapper = create_instance(self.settings["solver_wrappers"][master_sol_index])
        self.other_solver_wrapper_list = []

        for index, sol_wrapper_param in enumerate(self.settings["solver_wrappers"]):
            if not index == master_sol_index:
                self.other_solver_wrapper_list.append(create_instance(sol_wrapper_param))


        # create mappers
        self.mapper_interface_input_list = []
        for _ in self.other_solver_wrapper_list:
            self.mapper_interface_input_list.append(create_instance(self.settings["mapper_interface_input"]))

        self.mapper_interface_output_list = []
        for _ in self.other_solver_wrapper_list:
            self.mapper_interface_output_list.append(create_instance(self.settings["mapper_interface_output"]))

        # run time
        self.run_time = 0.0

    def initialize(self):
        super().initialize()

        self.master_solver_wrapper.initialize()
        for sol_wrapper in self.other_solver_wrapper_list:
            sol_wrapper.initialize()

        # create input mapper
        for index, other_sol_wrapper in enumerate(self.other_solver_wrapper_list):
            master_interface_input = self.master_solver_wrapper.get_interface_input().copy()
            other_interface_input = other_sol_wrapper.get_interface_input().copy()
            mapper_interface_input = self.mapper_interface_input_list[index]
            mapper_interface_input.initialize(master_interface_input, other_interface_input)

        # create output mapper
        for index, sol_wrapper in enumerate(self.other_solver_wrapper_list):
            interface_output_master = self.master_solver_wrapper.get_interface_output().copy()
            interface_output_other = sol_wrapper.get_interface_output().copy()
            mapper_interface_output = self.mapper_interface_output_list[index]
            mapper_interface_output.initialize(interface_output_other, interface_output_master)


    def initialize_solution_step(self):
        super().initialize_solution_step()
        self.master_solver_wrapper.initialize_solution_step()
        for other_sol_wrapper in self.other_solver_wrapper_list:
            other_sol_wrapper.initialize_solution_step()

    @tools.time_solve_solution_step
    def solve_solution_step(self, master_interface_input):
        self.master_interface_output = self.master_solver_wrapper.solve_solution_step(master_interface_input)
        for index, other_sol_wrapper in enumerate(self.other_solver_wrapper_list):
            other_interface_input = other_sol_wrapper.get_interface_input()
            mapper_interface_input = self.mapper_interface_input_list[index]
            mapper_interface_input(master_interface_input.copy(), other_interface_input)
            other_interface_output = other_sol_wrapper.solve_solution_step(other_interface_input)
            mapper_interface_output = self.mapper_interface_output_list[index]
            temp_master_interface_output = self.master_interface_output.copy()
            mapper_interface_output(other_interface_output, temp_master_interface_output)
            self.master_interface_output = self.master_interface_output + temp_master_interface_output

        return self.master_interface_output

    def finalize_solution_step(self):
        super().finalize_solution_step()

        self.master_solver_wrapper.finalize_solution_step()
        for other_sol_wrapper in self.other_solver_wrapper_list:
            other_sol_wrapper.finalize_solution_step()

    def finalize(self):
        super().finalize()

        self.master_solver_wrapper.finalize()
        for other_sol_wrapper in self.other_solver_wrapper_list:
            other_sol_wrapper.finalize()

        for mapper_interface_input in self.mapper_interface_input_list:
            mapper_interface_input.finalize()

        for mapper_interface_output in self.mapper_interface_output_list:
            mapper_interface_output.finalize()


    def output_solution_step(self):
        super().output_solution_step()

        for sol_wrapper in self.solver_wrapper_list:
            sol_wrapper.output_solution_step()

        for mapper_interface_input in self.mapper_interface_input_list:
            mapper_interface_input.output_solution_step()

        for mapper_interface_output in self.mapper_interface_output_list:
            mapper_interface_output.output_solution_step()

    def get_interface_input(self):
        # does not contain most recent data
        return self.master_solver_wrapper.get_interface_input()

    def get_interface_output(self):

        self.master_interface_output = self.master_solver_wrapper.get_interface_output()
        for index, other_sol_wrapper in self.other_solver_wrapper_list:
            other_interface_output = other_sol_wrapper.get_interface_output()
            mapper_interface_output = self.mapper_interface_output_list[index]
            temp_master_interface_output = self.master_interface_output.copy()
            mapper_interface_output(other_interface_output, temp_master_interface_output)
            self.master_interface_output += temp_master_interface_output

        return self.master_interface_output

    def print_components_info(self, pre):
        pass
