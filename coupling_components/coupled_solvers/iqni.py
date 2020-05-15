from coconut.coupling_components.tools import CreateInstance
from coconut.coupling_components.coupled_solvers.gauss_seidel import CoupledSolverGaussSeidel


def Create(parameters):
    return CoupledSolverIQNI(parameters)


class CoupledSolverIQNI(CoupledSolverGaussSeidel):
    def __init__(self, parameters):
        super().__init__(parameters)

        self.model = CreateInstance(self.parameters["settings"]["model"])
        self.omega = self.settings["omega"].GetDouble()

    def Initialize(self):
        super().Initialize()

        self.model.size_in = self.model.size_out = self.x.GetNumpyArray().shape[0]
        self.model.out = self.x.deepcopy()
        self.model.Initialize()
        self.components += [self.model]

    def SolveSolutionStep(self):
        # Initial value
        self.x = self.predictor.Predict(self.x)
        # First coupling iteration
        self.y = self.solver_wrappers[0].SolveSolutionStep(self.x)
        xt = self.solver_wrappers[1].SolveSolutionStep(self.y)
        r = xt - self.x
        self.model.Add(r, xt)
        self.FinalizeIteration(r)
        # Coupling iteration loop
        while not self.convergence_criterion.IsSatisfied():
            if not self.model.IsReady():
                dx = self.omega * r
            else:
                dr = -1 * r
                dx = self.model.Predict(dr) - dr
            self.x += dx
            self.y = self.solver_wrappers[0].SolveSolutionStep(self.x)
            xt = self.solver_wrappers[1].SolveSolutionStep(self.y)
            r = xt - self.x
            self.model.Add(r, xt)
            self.FinalizeIteration(r)

