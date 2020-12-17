from coconut import data_structure
from coconut.coupling_components import tools
from coconut.coupling_components.component import Component
from coconut.coupling_components.coupled_solvers.gauss_seidel import CoupledSolverGaussSeidel

import time
import os
import numpy as np


def create(parameters):
    return CoupledSolverTestSingleSolver(parameters)


class CoupledSolverTestSingleSolver(CoupledSolverGaussSeidel):

    def __init__(self, parameters):
        """"Should only initialize the solver that is to be tested"""
        Component.__init__(self)

        self.parameters = parameters
        self.settings = parameters.get("settings", None)  # settings is optional as long as the necessary parameters
        # are in test_settings

        if "test_settings" not in self.parameters.keys():  # requires a new parameter input "test_settings"
            raise KeyError('The coupled_solver "test_single_solver" requires "test_settings" which was not detected.')
        self.test_settings = parameters["test_settings"]
        self.solver_index = self.test_settings["solver_index"]  # solver to be tested; starts at 0
        self.test_class = self.test_settings.get("test_class", None)

        # copy value of settings to test_settings (test_settings are prioritized)
        if self.settings is not None:
            updated_settings = self.settings.update(self.test_settings)
            self.test_settings = updated_settings.copy()

        # delta_t and timestep_start
        self.timestep_start = self.test_settings.get("timestep_start", 0)
        self.test_settings.AddEmptyValue("timestep_start")
        self.test_settings.SetInt("timestep_start", self.timestep_start)
        self.n = self.timestep_start
        self.delta_t = self.test_settings["delta_t"]
        tools.print_info(f"Using delta_t = {self.delta_t} and timestep_start = {self.timestep_start}")

        self.predictor = DummyComponent()
        self.convergence_criterion = DummyComponent()
        # solver wrapper settings
        parameters = self.parameters["solver_wrappers"][self.solver_index]
        if parameters["type"] == "solver_wrappers.mapped":
            parameters = parameters["settings"]["solver_wrapper"]  # for mapped solver: the solver_wrapper itself tested
        settings = parameters["settings"]

        orig_wd = settings["working_directory"]  # working directory changed to a test_working_directory
        i = 0
        while os.path.exists(f"{orig_wd}_test{i}"):
            i += 1
        cur_wd = f"{orig_wd}_test{i}"
        settings.SetString("working_directory", cur_wd)
        os.system(f"cp -r {orig_wd} {cur_wd}")
        tools.print_info(f"{cur_wd} is the working_directory for the test\nCopying {orig_wd} to {cur_wd} \n")

        for key in ["timestep_start", "delta_t"]:  # add delta_t and timestep_start to solver_wrapper settings
            if key in settings:
                tools.print_info(f'WARNING: parameter "{key}" is defined multiple times in JSON file', layout='warning')
            settings[key] = self.test_settings[key]

        self.solver_wrapper = tools.create_instance(parameters)
        self.solver_wrappers = [self.solver_wrapper]  # used for printing summary

        self.components = [self.solver_wrapper]  # will only contain 1 solver wrapper

        # initialize test_class
        interface_input = self.solver_wrapper.interface_input
        if self.test_class is None:
            self.dummy_solver = None
            tools.print_info("No test class specified, zero input will be used")
            for model_part_name, variable in interface_input.model_part_variable_pairs:
                if data_structure.variables_dimensions[variable] == 1:
                    tools.print_info(f"\t0 is used as {variable} input to {model_part_name}")
                elif data_structure.variables_dimensions[variable] == 3:
                    tools.print_info(f"\t[0 0 0] is used as {variable} input to {model_part_name}")
        else:
            if not os.path.isfile('dummy_solver.py'):
                raise ModuleNotFoundError(f"Test class specified, but no file named dummy_solver.py in {os.getcwd()}")
            module = __import__('dummy_solver')
            if not hasattr(module, self.test_class):
                raise NameError(f"Module dummy_solver has no class {self.test_class}")
            self.dummy_solver = getattr(module, self.test_class)()
            tools.print_info(f"The functions from {self.test_class} will be used to calculate the following inputs:")
            for model_part_name, variable in interface_input.model_part_variable_pairs:
                if data_structure.variables_dimensions[variable] == 1:
                    tools.print_info(f"\t{variable} [Scalar] on {model_part_name}")
                elif data_structure.variables_dimensions[variable] == 3:
                    tools.print_info(f"\t{variable} [3D array] on {model_part_name}")
        tools.print_info()

        self.x = None
        self.y = None
        self.iteration = None  # iteration
        self.solver_level = 0  # 0 is main solver (time step is printed)

        self.start_time = None
        self.elapsed_time = None
        self.iterations = []
        self.save_results = self.test_settings.get("save_results", False)
        if self.save_results:
            self.complete_solution_x = None
            self.complete_solution_y = None
            self.residual = []
            self.case_name = self.test_settings.get("name", "results")  # case name
            self.case_name += "_" + cur_wd

    def initialize(self):
        Component.initialize(self)

        self.solver_wrapper.initialize()

        # initialize variables
        if self.solver_index == 1:
            self.x = self.solver_wrapper.get_interface_output()
            self.y = self.solver_wrapper.get_interface_input()
        else:
            self.x = self.solver_wrapper.get_interface_input()
            self.y = self.solver_wrapper.get_interface_output()

        if self.save_results:
            self.complete_solution_x = self.x.get_interface_data().reshape(-1, 1)
            self.complete_solution_y = self.y.get_interface_data().reshape(-1, 1)
        self.start_time = time.time()

    def solve_solution_step(self):
        interface_input = self.solver_wrapper.interface_input
        # generation of the input data
        if self.dummy_solver is not None:
            for model_part_name, variable in interface_input.model_parts_variables:
                model_part = interface_input.get_model_part(model_part_name)
                data = [getattr(self.dummy_solver, f"calculate_{variable}")(model_part.X0[i], model_part.Y0[i],
                                                                            model_part.Z0[i], self.n)
                        for i in range(model_part.size)]
                interface_input.set_varialbe_data(model_part_name, variable, np.array(data))
        # store data in self.x and self.y
        if self.solver_index == 1:
            self.y = interface_input
            self.x = self.solver_wrapper.solve_solution_step(interface_input)
        else:
            self.x = interface_input
            self.y = self.solver_wrapper.solve_solution_step(interface_input)
        self.finalize_iteration(self.x * 0)

    def print_header(self):
        if self.n == self.timestep_start + 1:
            header = f"════════════════════════════════════════════════════════════════════════════════\n" \
                f"{'Time step':<16}{'Norm x':<28}{'Norm y':<28}"
            tools.print_info(header, flush=True)

    def print_iteration_info(self, r):
        info = f"{self.n:<16d}{self.x.norm():<28.17e}{self.y.norm():<28.17e}"
        tools.print_info(' │' * self.solver_level, info, flush=True)


class DummyComponent:
    def update(self, x):
        pass
