from coconut.coupling_components.coupled_solvers.coupled_solver import CoupledSolver


def create(parameters):
    return CoupledSolverGaussSeidel(parameters)


class CoupledSolverGaussSeidel(CoupledSolver):

    def solve_solution_step(self):
        # initial value
        self.x = self.predictor.predict(self.x)
        print("Ini x: ", self.x)
        # first coupling iteration
        self.y = self.solver_wrappers[0].solve_solution_step(self.x.copy()).copy()
        print("y: ", self.y)
        xt = self.solver_wrappers[1].solve_solution_step(self.y.copy()).copy()
        r = xt - self.x
        self.finalize_iteration(r)
        # coupling iteration loop
        while not self.convergence_criterion.is_satisfied():
            self.x += r
            print("Updated x: ", self.x)
            self.y = self.solver_wrappers[0].solve_solution_step(self.x.copy()).copy()
            print("Updated y: ", self.y)
            xt = self.solver_wrappers[1].solve_solution_step(self.y.copy()).copy()
            r = xt - self.x
            self.finalize_iteration(r)
