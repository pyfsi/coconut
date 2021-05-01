from coconut.coupling_components.coupled_solvers.gauss_seidel import CoupledSolverGaussSeidel

import numpy as np


def create(parameters):
    return CoupledSolverAITKEN(parameters)


class CoupledSolverAITKEN(CoupledSolverGaussSeidel):
    def __init__(self, parameters):
        super().__init__(parameters)

        self.settings = parameters["settings"]
        self.omega_max = self.settings["omega_max"]

        self.omega = self.omega_max
        self.added = False
        self.rcurr = None

    def initialize(self):
        super().initialize()

        if self.restart:  # restart
            self.omega = self.restart_data['omega']

    def initialize_solution_step(self):
        super().initialize_solution_step()

        self.added = False
        self.rcurr = None

    def predict(self, r_in):
        r = r_in.get_interface_data()
        # calculate return value if sufficient data available
        if not self.added:
            raise RuntimeError("No information to predict")
        dx = self.omega * r
        dx_out = r_in.copy()
        dx_out.set_interface_data(dx)
        return dx_out

    def update(self, x_in, xt_in):
        x = x_in.get_interface_data().reshape(-1, 1)
        xt = xt_in.get_interface_data().reshape(-1, 1)
        r = xt - x
        rprev = self.rcurr
        self.rcurr = r
        if self.added:
            # Aitken relaxation
            # update omega
            self.omega *= -float(rprev.T @ (r - rprev) / np.linalg.norm(r - rprev, 2) ** 2)
        else:
            # set first value of omega in a time step
            self.omega = np.sign(self.omega) * min(abs(self.omega), self.omega_max)
            self.added = True

    def solve_solution_step(self):
        # initial value
        self.x = self.predictor.predict(self.x)
        # first coupling iteration
        self.y = self.solver_wrappers[0].solve_solution_step(self.x)
        xt = self.solver_wrappers[1].solve_solution_step(self.y)
        r = xt - self.x
        self.update(self.x, xt)
        self.finalize_iteration(r)
        # coupling iteration loop
        while not self.convergence_criterion.is_satisfied():
            self.x += self.predict(r)
            self.y = self.solver_wrappers[0].solve_solution_step(self.x)
            xt = self.solver_wrappers[1].solve_solution_step(self.y)
            r = xt - self.x
            self.update(self.x, xt)
            self.finalize_iteration(r)

    def is_ready(self):
        return self.added

    def add_restart_data(self, restart_data):
        return restart_data.update({'omega': self.omega})
