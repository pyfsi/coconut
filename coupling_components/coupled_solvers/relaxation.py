from coconut.coupling_components.coupled_solvers.gauss_seidel import CoupledSolverGaussSeidel


def create(parameters):
    return CoupledSolverRelaxation(parameters)


class CoupledSolverRelaxation(CoupledSolverGaussSeidel):
    def __init__(self, parameters):
        super().__init__(parameters)

        self.omega = self.settings["omega"]

    def solve_solution_step(self):
        # initial value
        self.x = self.predictor.predict(self.x)
        # first coupling iteration
        self.y = self.solver_wrappers[0].solve_solution_step(self.x)
        xt = self.solver_wrappers[1].solve_solution_step(self.y)
        r = xt - self.x
        self.finalize_iteration(r)
        # coupling iteration loop
        while not self.convergence_criterion.is_satisfied():
            self.x += self.omega * r
            self.y = self.solver_wrappers[0].solve_solution_step(self.x)
            xt = self.solver_wrappers[1].solve_solution_step(self.y)
            r = xt - self.x
            self.finalize_iteration(r)
