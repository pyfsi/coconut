from coconut import data_structure
from coconut.coupling_components import tools
from coconut.coupling_components.component import Component
from coconut.coupling_components.coupled_solvers.gauss_seidel import CoupledSolverGaussSeidel

import time
import os
import numpy as np


def Create(parameters):
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
        tools.Print(f"Using delta_t = {self.delta_t} and timestep_start = {self.timestep_start}")

        self.predictor = DummyComponent()
        self.convergence_criterion = DummyComponent()
        # solver wrapper settings
        parameters = self.parameters["solver_wrappers"][self.solver_index]
        if parameters["type"].GetString() == "solver_wrappers.mapped":
            parameters = parameters["settings"]["solver_wrapper"]  # for mapped solver: the solver_wrapper itself tested
        settings = parameters["settings"]

        orig_wd = settings["working_directory"].GetString()  # working directory changed to a test_working_directory
        i = 0
        while os.path.exists(f"{orig_wd}_test{i}"):
            i += 1
        cur_wd = f"{orig_wd}_test{i}"
        settings.SetString("working_directory", cur_wd)
        os.system(f"cp -r {orig_wd} {cur_wd}")
        tools.Print(f"{cur_wd} is the working_directory for the test\nCopying {orig_wd} to {cur_wd} \n")

        for key in ["timestep_start", "delta_t"]:  # add delta_t and timestep_start to solver_wrapper settings
            if settings.Has(key):
                tools.Print(f'WARNING: parameter "{key}" is defined multiple times in JSON file', layout='warning')
                settings.RemoveValue(key)
            settings.AddValue(key, self.test_settings[key])

        self.solver_wrapper = tools.CreateInstance(parameters)
        self.solver_wrappers = [self.solver_wrapper]  # used for printing summary

        self.components = [self.solver_wrapper]  # will only contain 1 solver wrapper

        # initialize test_class
        interface_input = self.solver_wrapper.interface_input
        if self.test_class is None:
            self.dummy_solver = None
            tools.Print("No test class specified, zero input will be used")
            for model_part_name, variable_names in interface_input.model_parts_variables:
                for variable_name in variable_names.list():
                    variable = vars(data_structure)[variable_name.GetString()]
                    if variable.Type() == "Double":
                        tools.Print(f"\t0 is used as {variable_name.GetString()} input to {model_part_name}")
                    elif variable.Type() == "Array":
                        tools.Print(f"\t[0 0 0] is used as {variable_name.GetString()} input to {model_part_name}")
        else:
            if not os.path.isfile('dummy_solver.py'):
                raise ModuleNotFoundError(f"Test class specified, but no file named dummy_solver.py in {os.getcwd()}")
            module = __import__('dummy_solver')
            if not hasattr(module, self.test_class):
                raise NameError(f"Module dummy_solver has no class {self.test_class}")
            self.dummy_solver = getattr(module, self.test_class)()
            tools.Print(f"The functions from {self.test_class} will be used to calculate the following inputs:")
            for model_part_name, variable_names in interface_input.model_parts_variables:
                for variable_name in variable_names.list():
                    variable = vars(data_structure)[variable_name.GetString()]
                    if variable.Type() == "Double":
                        tools.Print(f"\t{variable_name.GetString()} [Scalar] on {model_part_name}")
                    elif variable.Type() == "Array":
                        tools.Print(f"\t{variable_name.GetString()} [3D array] on {model_part_name}")
        tools.Print()

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

    def Initialize(self):
        Component.Initialize(self)

        self.solver_wrapper.Initialize()

        # Initialize variables
        if self.solver_index == 1:
            self.x = self.solver_wrapper.GetInterfaceOutput()
            self.y = self.solver_wrapper.GetInterfaceInput()
        else:
            self.x = self.solver_wrapper.GetInterfaceInput()
            self.y = self.solver_wrapper.GetInterfaceOutput()

        if self.save_results:
            self.complete_solution_x = self.x.GetNumpyArray().reshape(-1, 1)
            self.complete_solution_y = self.y.GetNumpyArray().reshape(-1, 1)
        self.start_time = time.time()

    def SolveSolutionStep(self):
        interface_input = self.solver_wrapper.interface_input
        # Generation of the input data
        if self.dummy_solver is not None:
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
            self.x = self.solver_wrapper.SolveSolutionStep(interface_input)
        else:
            self.x = interface_input
            self.y = self.solver_wrapper.SolveSolutionStep(interface_input)
        self.FinalizeIteration(self.x * 0)

    def PrintHeader(self):
        if self.n == self.timestep_start + 1:
            header = f"════════════════════════════════════════════════════════════════════════════════\n" \
                f"{'Time step':<16}{'Norm x':<28}{'Norm y':<28}"
            tools.Print(header)

    def PrintInfo(self, r):
        normx = np.linalg.norm(self.x.GetNumpyArray())
        normy = np.linalg.norm(self.y.GetNumpyArray())
        info = f"{self.n:<16d}{normx:<28.17e}{normy:<28.17e}"
        tools.Print(' │' * self.solver_level, info)


class DummyComponent:
    def Update(self, x):
        pass
