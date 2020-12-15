from coconut import data_structure
from coconut.coupling_components import tools
from coconut.coupling_components.component import Component
from coconut.coupling_components.coupled_solvers.gauss_seidel import CoupledSolverGaussSeidel

import time
import os
import sys
import numpy as np

try:
    from dummy_solver import *
    dummy_solver_present = True
except ImportError:
    dummy_solver_present = False


def create(parameters):
    return CoupledSolverTestSingleSolver(parameters)


class CoupledSolverTestSingleSolver(CoupledSolverGaussSeidel):

    def __init__(self, parameters):
        """"Should only initialize the solver that is to be tested"""
        Component.__init__(self)

        self.parameters = parameters
        self.settings = parameters["settings"] if self.parameters.Has("settings") else None  # settings is optional as long as the necessary parameters are in test_settings

        if "test_settings" not in self.parameters.keys():  # requires a new parameter input "test_settings"
            raise KeyError('The coupled_solver "test_single_solver" requires "test_settings" which was not detected.')
        self.test_settings = parameters["test_settings"]
        self.solver_index = self.test_settings["solver_index"].GetInt()  # solver to be tested; starts at 0
        self.test_class = self.test_settings["test_class"].GetString() if self.test_settings.Has("test_class") else None
        if self.test_class == "None":
            self.test_class = None

        # copy value of settings to test_settings (test_settings are prioritized)
        if self.settings is not None:
            self.test_settings.AddMissingParameters(self.settings)

        # delta_t and timestep_start
        self.timestep_start = self.test_settings["timestep_start"].GetInt() if self.test_settings.Has("timestep_start")\
            else 0
        self.test_settings.AddEmptyValue("timestep_start")
        self.test_settings.SetInt("timestep_start", self.timestep_start)
        self.n = self.timestep_start
        self.delta_t = self.test_settings["delta_t"].GetDouble()
        tools.print_info(f"Using delta_t = {self.delta_t} and timestep_start = {self.timestep_start}")

        self.predictor = DummyComponent()
        self.convergence_criterion = DummyComponent()
        # solver wrapper settings
        parameters = self.parameters["solver_wrappers"][self.solver_index]
        if parameters["type"].GetString() == "solver_wrappers.mapped":
            parameters = parameters["settings"]["solver_wrapper"]  # for mapped solver: the solver_wrapper itself tested
        settings = parameters["settings"]

        for key in ["timestep_start", "delta_t"]:  # add delta_t and timestep_start to solver_wrapper settings
            if settings.Has(key):
                tools.print_info(f'WARNING: parameter "{key}" is defined multiple times in JSON file', layout='warning')
                settings.RemoveValue(key)
            settings.AddValue(key, self.test_settings[key])

        self.solver_wrapper = tools.CreateInstance(parameters)
        self.solver_wrappers = [self.solver_wrapper]  # used for printing summary

        # working directory will be changed to a test_working_directory
        orig_wd = settings["working_directory"].GetString()
        i = 0
        while os.path.exists(f"{orig_wd}_test{i}"):
            i += 1
        cur_wd = f"{orig_wd}_test{i}"
        settings.SetString("working_directory", cur_wd)
        os.system(f"cp -r {orig_wd} {cur_wd}")
        tools.print_info(f"{cur_wd} is the working_directory for the test\nCopying {orig_wd} to {cur_wd} \n")

        self.components = [self.solver_wrapper]  # will only contain 1 solver wrapper

        interface_input = self.solver_wrapper.interface_input
        if self.test_class is None:
            self.zero_input = True
            tools.print_info("No test class specified, zero input will be used")
            for model_part_name, variable_names in interface_input.model_parts_variables:
                for variable_name in variable_names.list():
                    variable = vars(data_structure)[variable_name.GetString()]
                    if variable.Type() == "Double":
                        tools.print_info(f"\t0 is used as {variable_name.GetString()} input to {model_part_name}")
                    elif variable.Type() == "Array":
                        tools.print_info(f"\t[0 0 0] is used as {variable_name.GetString()} input to {model_part_name}")
        elif not dummy_solver_present:
            raise FileNotFoundError(f"Test class specified, but no file named dummy_solver.py in {os.getcwd()}")
        else:
            self.zero_input = False
            tools.print_info(f"The functions from {self.test_class} will be used to calculate the following inputs:")
            self.dummy_solver = getattr(sys.modules[__name__], self.test_class)()
            for model_part_name, variable_names in interface_input.model_parts_variables:
                for variable_name in variable_names.list():
                    variable = vars(data_structure)[variable_name.GetString()]
                    if variable.Type() == "Double":
                        tools.print_info(f"\t{variable_name.GetString()} [Scalar] on {model_part_name}")
                    elif variable.Type() == "Array":
                        tools.print_info(f"\t{variable_name.GetString()} [3D array] on {model_part_name}")
        tools.print_info()

        self.x = None
        self.y = None
        self.iteration = None  # Iteration
        self.solver_level = 0  # 0 is main solver (time step is printed)

        self.start_time = None
        self.elapsed_time = None
        self.iterations = []
        self.save_results = self.test_settings["save_results"].GetBool() if self.test_settings.Has("save_results") \
            else False
        if self.save_results:
            self.complete_solution_x = None
            self.complete_solution_y = None
            self.residual = []
            self.case_name = self.test_settings["name"].GetString() if self.test_settings.Has("name") else "results"  # Case name
            self.case_name += "_" + cur_wd

    def initialize(self):
        Component.initialize(self)

        self.solver_wrapper.initialize()

        # Initialize variables
        if self.solver_index == 1:
            self.x = self.solver_wrapper.get_interface_output()
            self.y = self.solver_wrapper.get_interface_input()
        else:
            self.x = self.solver_wrapper.get_interface_input()
            self.y = self.solver_wrapper.get_interface_output()

        if self.save_results:
            self.complete_solution_x = self.x.GetNumpyArray().reshape(-1, 1)
            self.complete_solution_y = self.y.GetNumpyArray().reshape(-1, 1)
        self.start_time = time.time()

    def solve_solution_step(self):
        interface_input = self.solver_wrapper.interface_input
        # Generation of the input data
        if not self.zero_input:
            for model_part_name, variable_names in interface_input.model_parts_variables:
                model_part = interface_input.model.GetModelPart(model_part_name)
                for variable_name in variable_names.list():
                    variable = vars(data_structure)[variable_name.GetString()]
                    for node in model_part.Nodes:
                        value = getattr(self.dummy_solver, f"calculate_{variable_name.GetString()}")(node.X0, node.Y0,
                                                                                                     node.Z0, self.n)
                        node.SetSolutionStepValue(variable, 0, value)
        # Store data in self.x and self.y
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
            tools.print_info(header)

    def print_iteration_info(self, r):
        info = f"{self.n:<16d}{self.x.norm():<28.17e}{self.y.norm():<28.17e}"
        tools.print_info(' │' * self.solver_level, info)


class DummyComponent:
    def update(self, x):
        pass
