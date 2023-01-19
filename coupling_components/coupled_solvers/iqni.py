from coconut.tools import create_instance, pass_on_parameters
from coconut.coupling_components.coupled_solvers.coupled_solver import CoupledSolver


def create(parameters):
    return CoupledSolverIQNI(parameters)


class CoupledSolverIQNI(CoupledSolver):
    def __init__(self, parameters):
        super().__init__(parameters)

        pass_on_parameters(self.settings, self.settings['model']['settings'], ('timestep_start', 'delta_t'))

        self.model = create_instance(self.parameters['settings']['model'])
        self.omega = self.settings['omega']

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
        self.y = self.solver_wrappers[0].solve_solution_step(self.x.copy()).copy()
        xt = self.solver_wrappers[1].solve_solution_step(self.y.copy()).copy()
        r = xt - self.x
        self.model.add(r.copy(), xt)
        self.finalize_iteration(r)
        # coupling iteration loop
        while not self.convergence_criterion.is_satisfied():
            if not self.model.is_ready():
                dx = self.omega * r
            else:
                dr = -1 * r
                dx = self.model.predict(dr.copy()) - dr
            self.x += dx
            self.y = self.solver_wrappers[0].solve_solution_step(self.x.copy()).copy()
            xt = self.solver_wrappers[1].solve_solution_step(self.y.copy()).copy()
            r = xt - self.x
            self.model.add(r.copy(), xt)
            self.finalize_iteration(r)

    def check_restart_data(self, restart_data):
        model_original = self.parameters['settings']['model']['type']
        model_new = restart_data['parameters']['settings']['model']['type']
        if model_original != model_new:
            raise ValueError(f'Restart not possible because model type changed:'
                             f'\n\toriginal: {model_original}\n\tnew: {model_new}')

    def add_restart_data(self, restart_data):
        return restart_data.update({'model': self.model})
