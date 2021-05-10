from coconut.tools import create_instance
from coconut.coupling_components.coupled_solvers.gauss_seidel import CoupledSolverGaussSeidel


def create(parameters):
    return CoupledSolverIQNI(parameters)


class CoupledSolverIQNI(CoupledSolverGaussSeidel):
    def __init__(self, parameters):
        super().__init__(parameters)

        self.model = create_instance(self.parameters["settings"]["model"])
        self.omega = self.settings["omega"]

    def initialize(self):
        super().initialize()

        if not self.restart:  # no restart
            self.model.size_in = self.model.size_out = self.x.size
            self.model.out = self.x.copy()
            self.model.initialize()
        else:  # restart
            self.model = self.restart_data['model']
        self.components += [self.model]

    def solve_solution_step(self):
        # initial value
        self.x = self.predictor.predict(self.x)
        # first coupling iteration
        self.y = self.solver_wrappers[0].solve_solution_step(self.x)
        xt = self.solver_wrappers[1].solve_solution_step(self.y)
        r = xt - self.x
        self.model.add(r, xt)
        self.finalize_iteration(r)
        # coupling iteration loop
        while not self.convergence_criterion.is_satisfied():
            if not self.model.is_ready():
                dx = self.omega * r
            else:
                dr = -1 * r
                dx = self.model.predict(dr) - dr
            self.x += dx
            self.y = self.solver_wrappers[0].solve_solution_step(self.x)
            xt = self.solver_wrappers[1].solve_solution_step(self.y)
            r = xt - self.x
            self.model.add(r, xt)
            self.finalize_iteration(r)

    def add_restart_check(self, restart_data):
        if self.parameters['settings']['model'] != restart_data['parameters']['settings']['model']:
            raise ValueError('Restart not possible because model changed')

    def add_restart_data(self, restart_data):
        return restart_data.update({'model': self.model})
