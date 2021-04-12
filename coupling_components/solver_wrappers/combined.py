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

        self.other_interface_output_list = []
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
        for index, sol_wrapper in enumerate(self.other_solver_wrapper_list):
            interface_input_from = self.master_solver_wrapper.get_interface_input().copy()
            interface_input_to = sol_wrapper.get_interface_input().copy()
            mapper_interface_input = self.mapper_interface_input_list[index]
            mapper_interface_input.initialize(interface_input_from, interface_input_to)

        # create output mapper
        for index, sol_wrapper in enumerate(self.other_solver_wrapper_list):
            interface_output_to = self.master_solver_wrapper.get_interface_output().copy()
            interface_output_from = sol_wrapper.get_interface_output().copy()
            mapper_interface_output = self.mapper_interface_output_list[index]
            mapper_interface_output.initialize(interface_output_from, interface_output_to)


    def initialize_solution_step(self):
        super().initialize_solution_step()

        for sol_wrapper in self.solver_wrapper_list:
            sol_wrapper.initialize_solution_step()

    @tools.time_solve_solution_step
    def solve_solution_step(self, master_interface_input):
        self.master_interface_output = self.master_solver_wrapper.solve_solution_step(master_interface_input)
        for index, sol_wrapper in self.other_solver_wrapper_list:
            other_interface_input = sol_wrapper.get_interface_input()
            mapper_interface_input = self.mapper_interface_input_list[index]
            mapper_interface_input(master_interface_input.copy(), other_interface_input)
            other_interface_output = sol_wrapper.solve_solution_step(other_interface_input)
            mapper_interface_output = self.mapper_interface_output_list[index]
            temp_master_interface_output = self.master_interface_output.copy()
            mapper_interface_output(other_interface_output, temp_master_interface_output)
            self.master_interface_output += temp_master_interface_output

        return self.master_interface_output

    def finalize_solution_step(self):
        super().finalize_solution_step()

        for sol_wrapper in self.solver_wrapper_list:
            sol_wrapper.finalize_solution_step()

    def finalize(self):
        super().finalize()

        for sol_wrapper in self.solver_wrapper_list:
            sol_wrapper.finalize()

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
        return self.master_interface_input

    def get_interface_output(self):
        interface_output_from = self.solver_wrapper.get_interface_output() #TODO: map the fields to master and send output
        self.mapper_interface_output(interface_output_from, self.interface_output_to)
        return self.interface_output_to

    def print_components_info(self, pre):
        tools.print_info(pre, "The component ", self.__class__.__name__, " maps the following solver wrapper:")
        pre = tools.update_pre(pre)
        self.solver_wrapper.print_components_info(pre + '├─')
        tools.print_info(pre, '├─', "Input mapper:")
        self.mapper_interface_input.print_components_info(pre + '│ └─')
        tools.print_info(pre, '└─', "Output mapper:")
        self.mapper_interface_output.print_components_info(pre + '  └─')
