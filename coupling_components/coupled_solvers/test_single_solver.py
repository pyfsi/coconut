from coconut import data_structure
from coconut.coupling_components import tools
from coconut.coupling_components.component import Component
from coconut.coupling_components.coupled_solvers.gauss_seidel import CoupledSolverGaussSeidel

import time
import os
import sys

try:
    from dummy_solver import *
    bool_test = True
except ImportError:
    bool_test = False
    print(f"No file named dummy_solver.py found, zero input will be used")


def Create(parameters):
    return CoupledSolverTestSingleSolver(parameters)


class CoupledSolverTestSingleSolver(CoupledSolverGaussSeidel):

    def __init__(self, parameters):
        """"Should only initialize the solver that is to be tested"""
        Component.__init__(self)

        self.parameters = parameters
        if "test_settings" not in self.parameters.keys():
            raise KeyError("The coupled_solver \"test_single_solver\" requires "
                           "\"test_settings\" which was not detected.")

        # settings is optional as long as the necessary parameters are in test_settings
        self.settings = parameters["settings"] if self.parameters.Has("settings") else None

        self.test_settings = parameters["test_settings"]  # requires a new parameter input "test_settings"
        self.solver_index = self.test_settings["solver_index"].GetInt()  # starts at 0
        self.test_class = self.test_settings["test_class"].GetString() if self.test_settings.Has("test_class") else None
        self.test_class = self.test_class if not self.test_class == "None" else None

        # copy value of settings to test_settings (test_settings are prioritized)
        if self.settings is not None:
            self.test_settings.AddMissingParameters(self.settings)

        # delta_t and timestep_start
        self.delta_t = self.test_settings["delta_t"].GetDouble()
        print(f"Using delta_t = {self.delta_t}")

        tools.Print("timestep_start set to 0 (Default for test_single_solver)")
        self.n = 0

        # add delta_t and timestep_start to solver_wrapper settings
        parameters = self.parameters["solver_wrappers"][self.solver_index]
        if parameters["type"].GetString() == "solver_wrappers.mapped":
            parameters = parameters["settings"]["solver_wrapper"]  # for mapped solver: the solver_wrapper itself tested
        settings = parameters["settings"]

        if settings.Has("delta_t"):
            tools.Print(f'WARNING: parameter "{"delta_t"}" is defined multiple times in JSON file', layout='warning')
            settings.RemoveValue("delta_t")
        settings.AddValue("delta_t", self.test_settings["delta_t"])

        settings.AddEmptyValue("timestep_start")
        settings.SetInt("timestep_start", 0)

        # working directory will be changed to a test_working_directory
        orig_wd = settings["working_directory"].GetString()
        i = 0
        while os.path.exists(f"{orig_wd}_test{i}"):
            i += 1
        cur_wd = f"{orig_wd}_test{i}"
        settings.SetString("working_directory", cur_wd)
        os.system(f"cp -r {orig_wd} {cur_wd}")
        print(f"{cur_wd} is the working_directory for the test\nCopying {orig_wd} to {cur_wd} \n")

        self.solver_wrapper = tools.CreateInstance(parameters)

        self.components = [self.solver_wrapper]  # will only contain 1 solver wrapper

        input = self.solver_wrapper.interface_input
        if self.test_class is None:
            print("No test class specified")
            for model_part_name, variable_names in input.model_parts_variables:
                for variable_name in variable_names.list():
                    variable = vars(data_structure)[variable_name.GetString()]
                    if variable.Type() is "Double":
                        print(f"\t0 is used as {variable_name.GetString()} input to {model_part_name}")
                    elif variable.Type() is "Array":
                        print(f"\t[0 0 0] is used as {variable_name.GetString()} input to {model_part_name}")
        else:
            print(f"The functions from the class {self.test_class} will be used to calculate the following inputs:")
            self.Test_object = getattr(sys.modules[__name__], self.test_class)
            for model_part_name, variable_names in input.model_parts_variables:
                for variable_name in variable_names.list():
                    variable = vars(data_structure)[variable_name.GetString()]
                    if variable.Type() is "Double":
                        print(f"\t{variable_name.GetString()} [Scalar] on {model_part_name}")
                    elif variable.Type() is "Array":
                        print(f"\t{variable_name.GetString()} [3D array] on {model_part_name}")
        print("")

        self.start_time = None
        self.stop_time = None
        self.iterations = [1]
        self.save_results = False  # Always False

    def Initialize(self):
        Component.Initialize(self)

        self.solver_wrapper.Initialize()
        # for component in self.components[1:]:
        #     component.Initialize()

        #
        # # Construct mappers if required
        # index_mapped = None
        # index_other = None
        # for i in range(2):
        #     type = self.parameters["solver_wrappers"][i]["type"].GetString()
        #     if type == "solver_wrappers.mapped":
        #         index_mapped = i
        #     else:
        #         index_other = i
        # # if index_other is None:
        # #     raise ValueError("Not both solvers may be mapped solvers.")
        # if index_mapped is not None:
        #     # Construct input mapper
        #     interface_input_from = self.solver_wrappers[index_other].GetInterfaceOutput()
        #     self.solver_wrappers[index_mapped].SetInterfaceInput(interface_input_from)
        #
        #     # Construct output mapper
        #     interface_output_to = self.solver_wrappers[index_other].GetInterfaceInput()
        #     self.solver_wrappers[index_mapped].SetInterfaceOutput(interface_output_to)

        # # Initialize variables
        # self.x = self.solver_wrappers[1].GetInterfaceOutput()
        # self.predictor.Initialize(self.x)
        #

        self.solver_level = 0  # 0 is main solver (time step is printed)

        self.start_time = time.time()

    def SolveSolutionStep(self):
        interface_input = self.solver_wrapper.interface_input

        # Generation of the input data
        if not self.test_class == "None" and bool_test:
            for model_part_name, variable_names in interface_input.model_parts_variables:
                model_part = interface_input.model.GetModelPart(model_part_name)
                for variable_name in variable_names.list():
                    variable = vars(data_structure)[variable_name.GetString()]
                    for node in model_part.Nodes:
                        value = getattr(self.Test_object, f"calculate_{variable_name.GetString()}")(self.Test_object,
                                                                                                    node.X0, node.Y0,
                                                                                                    node.Z0, self.n)
                        node.SetSolutionStepValue(variable, 0, value)

        interface_output = self.solver_wrapper.SolveSolutionStep(interface_input)

    def FinalizeSolutionStep(self):
        Component.FinalizeSolutionStep(self)

        for component in self.components:
            component.FinalizeSolutionStep()
