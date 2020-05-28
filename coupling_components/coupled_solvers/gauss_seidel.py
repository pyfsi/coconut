from coconut.coupling_components import tools
from coconut.coupling_components.tools import CreateInstance
from coconut.coupling_components.component import Component

import numpy as np
import time
import pickle


def Create(parameters):
    return CoupledSolverGaussSeidel(parameters)


class CoupledSolverGaussSeidel(Component):
    def __init__(self, parameters):
        super().__init__()

        self.parameters = parameters
        self.settings = parameters["settings"]

        self.timestep_start = self.settings["timestep_start"].GetInt()  # Time step where calculation is started
        self.n = self.timestep_start  # Time step
        self.delta_t = self.settings["delta_t"].GetDouble()  # Time step size

        self.predictor = CreateInstance(self.parameters["predictor"])
        self.convergence_criterion = CreateInstance(self.parameters["convergence_criterion"])
        self.solver_wrappers = []
        for index in range(2):
            # Add timestep_start and delta_t to solver_wrapper settings
            parameters = self.parameters["solver_wrappers"][index]
            if parameters["type"].GetString() == "solver_wrappers.mapped":
                settings = parameters["settings"]["solver_wrapper"]["settings"]
            else:
                settings = parameters["settings"]

            for key in ["timestep_start", "delta_t"]:
                if settings.Has(key):
                    tools.Print(f'WARNING: parameter "{key}" is defined multiple times in JSON file', layout='warning')
                    settings.RemoveValue(key)
                settings.AddValue(key, self.settings[key])

            self.solver_wrappers.append(CreateInstance(parameters))

        self.components = [self.predictor, self.convergence_criterion, self.solver_wrappers[0], self.solver_wrappers[1]]

        self.x = None
        self.y = None
        self.iteration = None  # Iteration

        self.start_time = None
        self.stop_time = None
        self.iterations = []
        if self.settings.Has("save_results"):
            # Set True in order to save for every iteration
            self.save_results = self.settings["save_results"].GetBool()
        else:
            self.save_results = False
        if self.save_results:
            self.complete_solution_x = None
            self.complete_solution_y = None
            self.residual = []
            if self.settings.Has("name"):
                self.case_name = self.settings["name"].GetString()  # case name
            else:
                self.case_name = "results"

    def Initialize(self):
        super().Initialize()

        for component in self.components[1:]:
            component.Initialize()

        # Construct mappers if required
        index_mapped = None
        index_other = None
        for i in range(2):
            type = self.parameters["solver_wrappers"][i]["type"].GetString()
            if type == "solver_wrappers.mapped":
                index_mapped = i
            else:
                index_other = i
        if index_other is None:
            raise ValueError("Not both solvers may be mapped solvers.")
        if index_mapped is not None:
            # Construct input mapper
            interface_input_from = self.solver_wrappers[index_other].GetInterfaceOutput()
            self.solver_wrappers[index_mapped].SetInterfaceInput(interface_input_from)

            # Construct output mapper
            interface_output_to = self.solver_wrappers[index_other].GetInterfaceInput()
            self.solver_wrappers[index_mapped].SetInterfaceOutput(interface_output_to)

        # Initialize variables
        self.x = self.solver_wrappers[1].GetInterfaceOutput()
        self.y = self.solver_wrappers[0].GetInterfaceOutput()
        self.predictor.Initialize(self.x)

        if self.save_results:
            self.complete_solution_x = self.x.GetNumpyArray().reshape(-1, 1)
            self.complete_solution_y = self.y.GetNumpyArray().reshape(-1, 1)
        self.start_time = time.time()

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()

        for component in self.components:
            component.InitializeSolutionStep()

        self.n += 1  # Increment time step
        self.iteration = 0

        # Print timestep
        out = f"=============================================" \
              f"===================================\n" \
              f"\tTime step {self.n}\n" \
              f"=============================================" \
              f"===================================\n" \
              f"Iteration\tNorm residual"
        tools.PrintInfo(out)

        if self.save_results:
            self.residual.append([])

    def SolveSolutionStep(self):
        # Initial value
        self.x = self.predictor.Predict(self.x)
        # First coupling iteration
        self.y = self.solver_wrappers[0].SolveSolutionStep(self.x)
        xt = self.solver_wrappers[1].SolveSolutionStep(self.y)
        r = xt - self.x
        self.FinalizeIteration(r)
        # Coupling iteration loop
        while not self.convergence_criterion.IsSatisfied():
            self.x += r
            self.y = self.solver_wrappers[0].SolveSolutionStep(self.x)
            xt = self.solver_wrappers[1].SolveSolutionStep(self.y)
            r = xt - self.x
            self.FinalizeIteration(r)

    def FinalizeIteration(self, r):
        self.iteration += 1
        self.convergence_criterion.Update(r)
        # Print iteration information
        norm = np.linalg.norm(r.GetNumpyArray())
        out = f"{self.iteration:<9d}\t{norm:<22.17e}"
        tools.PrintInfo(out)

        if self.save_results:
            self.residual[self.n - 1].append(np.linalg.norm(r.GetNumpyArray()))

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()

        self.iterations.append(self.iteration)
        if self.save_results:
            self.complete_solution_x = np.hstack((self.complete_solution_x, self.x.GetNumpyArray().reshape(-1, 1)))
            self.complete_solution_y = np.hstack((self.complete_solution_y, self.y.GetNumpyArray().reshape(-1, 1)))

        self.predictor.Update(self.x)
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

        self.stop_time = time.time()
        elapsed_time = self.stop_time - self.start_time
        print(f"\nElapsed time: {elapsed_time:0.3f}s\n"
              f"Average number of iterations per time step: {np.array(self.iterations).mean():0.2f}")
        if self.save_results:
            output = {"solution_x": self.complete_solution_x, "solution_y": self.complete_solution_y,
                      "interface_x": self.x, "interface_y": self.y, "iterations": self.iterations,
                      "time": elapsed_time, "residual": self.residual, "delta_t": self.delta_t,
                      "timestep_start": self.timestep_start}
            pickle.dump(output, open(self.case_name, 'wb'))

    def Check(self):
        super().Check()

        for component in self.components:
            component.Check()

    def PrintInfo(self, indent):
        tools.Print('\n', '\t' * indent, "The coupled solver ", self.__class__.__name__,
                    " has the following components:")
        for component in self.components:
            component.PrintInfo(indent + 1)
