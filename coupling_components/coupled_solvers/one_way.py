from coconut.tools import print_info
from coconut.coupling_components.coupled_solvers.coupled_solver import CoupledSolver
import numpy as np

def create(parameters):
    return CoupledSolverOneWay(parameters)


class CoupledSolverOneWay(CoupledSolver):
    def __init__(self, parameters):
        super().__init__(parameters)
        print_info('CoupledSolverOneWay is chosen: convergence criterion and predictor are ignored', layout='info')

    def solve_solution_step(self):
        self.x = self.predictor.predict(self.x)
        self.y = self.solver_wrappers[0].solve_solution_step(0*self.x.copy()).copy() # set interface to zero: no effect of other solver
        xt = self.solver_wrappers[1].solve_solution_step(self.y.copy())
        r = xt - self.x
        self.x = xt  # for storing resulting interface
        self.finalize_iteration(r)
