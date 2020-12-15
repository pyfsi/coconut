from coconut.coupling_components.tools import CreateInstance
from coconut.coupling_components.coupled_solvers.gauss_seidel import CoupledSolverGaussSeidel


def create(parameters):
    return CoupledSolverIQNI(parameters)


class CoupledSolverIQNI(CoupledSolverGaussSeidel):
    def __init__(self, parameters):
        super().__init__(parameters)

        self.model = CreateInstance(self.parameters["settings"]["model"])
        self.omega = self.settings["omega"].GetDouble()

    def initialize(self):
        super().initialize()

        self.model.size_in = self.model.size_out = self.x.GetNumpyArray().shape[0]
        self.model.out = self.x.deepcopy()
        self.model.initialize()
        self.components += [self.model]

    def solve_solution_step(self):
        # Initial value
        self.x = self.predictor.Predict(self.x)
        # First coupling iteration
        self.y = self.solver_wrappers[0].solve_solution_step(self.x)
        xt = self.solver_wrappers[1].solve_solution_step(self.y)
        r = xt - self.x
        self.model.Add(r, xt)
        self.finalize_iteration(r)
        # Coupling iteration loop
        while not self.convergence_criterion.is_satisfied():
            if not self.model.IsReady():
                dx = self.omega * r
            else:
                dr = -1 * r
                dx = self.model.Predict(dr) - dr
            self.x += dx
            self.y = self.solver_wrappers[0].solve_solution_step(self.x)
            xt = self.solver_wrappers[1].solve_solution_step(self.y)
            r = xt - self.x
            self.model.Add(r, xt)
            self.finalize_iteration(r)
