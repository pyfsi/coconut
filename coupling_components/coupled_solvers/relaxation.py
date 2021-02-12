from coconut.coupling_components.coupled_solvers.gauss_seidel import CoupledSolverGaussSeidel


def Create(parameters):
    return CoupledSolverRelaxation(parameters)


class CoupledSolverRelaxation(CoupledSolverGaussSeidel):
    def __init__(self, parameters):
        super().__init__(parameters)

        self.omega = self.settings["omega"].GetDouble()

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
            self.x += self.omega * r
            self.y = self.solver_wrappers[0].SolveSolutionStep(self.x)
            xt = self.solver_wrappers[1].SolveSolutionStep(self.y)
            r = xt - self.x
            self.FinalizeIteration(r)
