from coconut import data_structure
from coconut.coupling_components import tools
from coconut.coupling_components.tools import CreateInstance
from coconut.coupling_components.component import Component

import numpy as np
import time
import os
import sys

try:
    from test_function import *
    bool_test = True
except ImportError:
    bool_test = False
    print(f"No file named test_function.py found, zero input will be used")

def Create(parameters):
    return CoupledSolverTestSingleSolver(parameters)

def Param_Priority(key, set1, set2): # Code will give priority to set1 but can also look in set2
    if set2:
        if key in set1.keys():
            return set1
        elif key in set2.keys():
            return set2
        else:
            raise KeyError(f"{key} is required but could not be found in set1 nor set2")
    else:
        if key in set1.keys():
            return set1
        else:
            raise KeyError(f"{key} is required but could not be found in set1 and there is no set2")

class CoupledSolverTestSingleSolver(Component):
    def __init__(self, parameters):
        '''Should only initialize the solver that is to be tested'''
        super().__init__()

        self.parameters = parameters
        if not ("test_settings" in self.parameters.keys()):
             raise KeyError("The coupled_solver \"test_single_solver\" requires \"test_settings\" which was not detected.")

        if ("settings" in self.parameters.keys()): #settings is optional as long as the necessary parameters are in test_settings
            self.settings = parameters["settings"]
        else:
            self.settings=None

        self.test_settings = parameters["test_settings"] #Requires a new parameter input "test_settings"
        self.solver_index = self.test_settings["solver_index"].GetInt() #Starts at 0

        self.test_class = self.test_settings["test_class"].GetString()



        ##delta_t
        self.delta_t = Param_Priority("delta_t",self.test_settings, self.settings)["delta_t"].GetDouble()
        print(f"Using delta_t = {self.delta_t}")

        # if "delta_t" in self.test_settings.keys():
        #     self.delta_t = self.test_settings["delta_t"].GetDouble()  # Time step size
        #     tools.Print(f"Using delta_t = {self.delta_t} from test_settings",layout="plain")
        # elif bool_settings and ("delta_t" in self.settings.keys()): #if no delta_t in test_settings use from settings
        #     self.delta_t = self.settings["delta_t"].GetDouble()  # Time step size
        #     tools.Print(f"Using delta_t = {self.delta_t} from settings",layout="plain")
        # else:
        #     raise KeyError("delta_t is required but could not be found in test_settings nor settings")

        #
        # self.n = self.settings["timestep_start"].GetInt()  # Time step
        # self.delta_t = self.settings["delta_t"].GetDouble()  # Time step size

        # self.predictor = CreateInstance(self.parameters["predictor"])
        # self.convergence_criterion = CreateInstance(self.parameters["convergence_criterion"])
        self.solver_wrappers = []
        tools.Print("timestep_start set to 0 (Default for test_single_solver)",layout="plain")
        self.n = 0


        # Add timestep_start and delta_t to solver_wrapper settings
        parameters = self.parameters["solver_wrappers"][self.solver_index]
        if parameters["type"].GetString() == "solver_wrappers.mapped":
            parameters = parameters["settings"]["solver_wrapper"] # for a mapped solver want to test the solver_wrapper itself
            settings = parameters["settings"]
        else:
            settings = parameters["settings"]

        for key in ["delta_t"]:
            if settings.Has(key):
                tools.Print(f'WARNING: parameter "{key}" is defined multiple times in JSON file', layout='warning')
                settings.RemoveValue(key)
            settings.AddEmptyValue(key)
            settings.SetDouble(key,self.delta_t)

        settings.AddEmptyValue("timestep_start")
        settings.SetInt("timestep_start",0)

        #Working directory will be changed to a test_working_directory
        orig_wd = settings["working_directory"].GetString()
        i=0
        while(os.path.exists(f"{orig_wd}_test{i}")):
                i+=1
        cur_wd = f"{orig_wd}_test{i}"
        settings.SetString("working_directory", cur_wd)
        os.system(f"cp -r {orig_wd} {cur_wd}")
        print(f"{cur_wd} is the working_directory for the test\nCopying {orig_wd} to {cur_wd} \n")



        self.solver_wrappers.append(CreateInstance(parameters))

        # self.components = [self.predictor, self.convergence_criterion, self.solver_wrappers[0], self.solver_wrappers[1]]
        self.components = [self.solver_wrappers[0]] #Will only contain 1 solver wrapper

        input = self.solver_wrappers[0].interface_input
        if not bool_test:
            print("No test_function.py was found:")
            for model_part_name, variable_names in input.model_parts_variables:
                # model_part = input.model.GetModelPart(model_part_name)
                for variable_name in variable_names.list():
                    variable = vars(data_structure)[variable_name.GetString()]
                    if variable.Type() is "Double":
                        print(f"\t0 is used as {variable_name.GetString()} input to {model_part_name}")
                    elif variable.Type() is "Array":
                        print(f"\t[0 0 0] is used as {variable_name.GetString()} input to {model_part_name}")

        elif self.test_class == "None":
            print("No test class specified")
            for model_part_name, variable_names in input.model_parts_variables:
                # model_part = input.model.GetModelPart(model_part_name)
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
                # model_part = input.model.GetModelPart(model_part_name)
                for variable_name in variable_names.list():
                    variable = vars(data_structure)[variable_name.GetString()]
                    if variable.Type() is "Double":
                        print(f"\t{variable_name.GetString()} [Scalar] on {model_part_name}")
                    elif variable.Type() is "Array":
                        print(f"\t{variable_name.GetString()} [3D array] on {model_part_name}")
            # print(input.model_parts_variables)
            # node_coords = input.GetInitialCoordinates()
            # print(node_coords.shape)
        print("")

        self.x = []
        self.iteration = None  # Iteration

        self.save_iterations = False  # Set True in order to save iteration related information
        if self.save_iterations:
            self.complete_solution = None
            self.complete_solution_y = None
            self.iterations = []
            self.start_time = None
            self.stop_time = None
            self.residual = []

    def Initialize(self):
        super().Initialize()

        self.solver_wrappers[0].Initialize()
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


        if self.save_iterations:
            self.start_time = time.time()

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()

        for component in self.components:
            component.InitializeSolutionStep()

        # print("solutionstep initialized")

        self.n += 1  # Increment time step
        self.iteration = 0

        # Print time step
        if not self.solver_level:
            out = f"════════════════════════════════════════════════════════════════════════════════\n" \
                  f"\tTime step {self.n}\n" \
                  f"════════════════════════════════════════════════════════════════════════════════\n" \
                  # f"Iteration\tNorm residual"
            tools.Print(out)

    def SolveSolutionStep(self):
        # # Initial value
        # self.x = self.predictor.Predict(self.x)
        # First coupling iteration
        # y = self.solver_wrappers[0].SolveSolutionStep(self.x)
        # P = self.solver_wrappers[0].interface_input.GetPythonList()
        input = self.solver_wrappers[0].interface_input


        #### Put manipulation of the input here by creating a data list ####
        if (not self.test_class == "None") and bool_test:
            for model_part_name, variable_names in input.model_parts_variables:
                model_part = input.model.GetModelPart(model_part_name)
                for variable_name in variable_names.list():
                    variable = vars(data_structure)[variable_name.GetString()]
                    data = []
                    for node in model_part.Nodes:
                        if variable.Type() is "Double":
                            value = getattr(self.Test_object,f"calculate_{variable_name.GetString()}")(node.X0,node.Y0,node.Z0,0)
                            node.SetSolutionStepValue(variable, 0, value)
                        elif variable.Type() is "Array":
                            value = getattr(self.Test_object,f"calculate_{variable_name.GetString()}")(node.X0,node.Y0,node.Z0,0)
                            node.SetSolutionStepValue(variable, 0, value)

        elif self.test_class == "None":
            print("No test class specified")
            for model_part_name, variable_names in input.model_parts_variables:
                # model_part = input.model.GetModelPart(model_part_name)
                for variable_name in variable_names.list():
                    variable = vars(data_structure)[variable_name.GetString()]
                    if variable.Type() is "Double":
                        print(f"\t0 is used as {variable_name.GetString()} input to {model_part_name}")
                    elif variable.Type() is "Array":
                        print(f"\t[0 0 0] is used as {variable_name.GetString()} input to {model_part_name}")
        else:
            print(f"The functions from the class {self.test_class} will be used to calculate the following inputs:")
            for model_part_name, variable_names in input.model_parts_variables:
                # model_part = input.model.GetModelPart(model_part_name)
                for variable_name in variable_names.list():
                    variable = vars(data_structure)[variable_name.GetString()]
                    if variable.Type() is "Double":
                        print(f"\t{variable_name.GetString()} [Scalar] on {model_part_name}")
                    elif variable.Type() is "Array":
                        print(f"\t{variable_name.GetString()} [3D array] on {model_part_name}")
            # print(input.model_parts_variables)
            # node_coords = input.GetInitialCoordinates()
            # print(node_coords.shape)
        ####################################################################


        output = self.solver_wrappers[0].SolveSolutionStep(input)
        r = 0
        self.FinalizeIteration(r)
        # Coupling iteration loop
        # while not self.convergence_criterion.IsSatisfied():
        #     self.x += r
        #     y = self.solver_wrappers[0].SolveSolutionStep(self.x)
        #     xt = self.solver_wrappers[1].SolveSolutionStep(y)
        #     r = xt - self.x
        #     self.FinalizeIteration(r)

    def FinalizeIteration(self, r):
        self.iteration += 1
        # self.convergence_criterion.Update(r)
        # # Print iteration information
        # norm = np.linalg.norm(r.GetNumpyArray())
        # out = f"{self.iteration:<9d}\t{norm:<22.17e}"
        # tools.PrintInfo(out)

        if self.save_iterations:
            self.residual[self.n - 1].append(np.linalg.norm(r.GetNumpyArray()))

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()

        # if self.save_iterations:
        #     timestep_solution = self.x.GetNumpyArray().reshape(-1, 1)
        #     y = self.solver_wrappers[0].SolveSolutionStep(self.x)
        #     timestep_solution_y = y.GetNumpyArray().reshape(-1, 1)
        #     if self.complete_solution is None:
        #         self.complete_solution = timestep_solution
        #         self.complete_solution_y = timestep_solution_y
        #     else:
        #         self.complete_solution = np.hstack((self.complete_solution, timestep_solution))
        #         self.complete_solution_y = np.hstack((self.complete_solution_y, timestep_solution_y))
        #     self.iterations.append(self.iteration)

        # self.predictor.Update(self.x)
        for component in self.components:
            component.FinalizeSolutionStep()

    def OutputSolutionStep(self):
        super().OutputSolutionStep()

        for component in self.components:
            component.OutputSolutionStep()

    def Finalize(self):
        super().Finalize()

        for component in self.components:
            component.Finalize()

        if self.save_iterations:
            self.stop_time = time.time()
            type = self.parameters["type"].GetString()
            if self.parameters["settings"].Has("model"):
                model = '_' + self.parameters["settings"]["model"]["type"].GetString()
                if self.parameters["settings"]["model"]["settings"].Has("q") and model == "_coupled_solvers.models.ls":
                    q = '_q' + str(self.parameters["settings"]["model"]["settings"]["q"].GetDouble())
                else:
                    q = ''
            else:
                model = ''
                q = ''
            if self.parameters["settings"].Has("surrogate"):
                sur = '_' + self.parameters["settings"]["surrogate"]["type"].GetString()[23:]
            else:
                sur = ''

            output_name = 'result.' + type + model + q + sur
            output = {"solution": self.complete_solution, "solution_y": self.complete_solution_y,
                      "iterations": self.iterations, "time": self.stop_time - self.start_time,
                      "residual": self.residual}
            np.save(output_name, output)

    def Check(self):
        super().Check()

        for component in self.components:
            component.Check()

    def PrintInfo(self, pre):
        tools.Print(pre, "The coupled solver ", self.__class__.__name__, " has the following components:")
        tools.PrintStructureInfo(pre, self.components)

