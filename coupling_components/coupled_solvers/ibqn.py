from coconut.tools import create_instance, pass_on_parameters
from coconut.coupling_components.coupled_solvers.coupled_solver import CoupledSolver


from scipy.sparse.linalg import gmres, LinearOperator


def create(parameters):
    return CoupledSolverIBQN(parameters)


class CoupledSolverIBQN(CoupledSolver):
    def __init__(self, parameters):
        super().__init__(parameters)

        for model in ('model_f', 'model_s'):
            pass_on_parameters(self.settings, self.settings[model]['settings'], ('timestep_start', 'delta_t'))

        self.model_f = create_instance(self.parameters['settings']['model_f'])
        self.model_s = create_instance(self.parameters['settings']['model_s'])
        self.omega = self.settings['omega']
        self.atol = self.settings['absolute_tolerance_gmres']
        self.rtol = self.settings['relative_tolerance_gmres']

        self.xtemp = self.ytemp = None
        self.dxtemp = self.dytemp = None
        self.u = self.w = None
        self.ready = None

    def initialize(self):
        super().initialize()

        self.dxtemp = self.x.copy()
        self.dytemp = self.y.copy()
        self.u = self.x.size
        self.w = self.y.size
        self.ready = False
        if not self.restart:  # no restart
            self.model_f.size_in = self.model_s.size_out = self.u
            self.model_f.size_out = self.model_s.size_in = self.w
            self.model_f.out = self.y.copy()
            self.model_s.out = self.x.copy()
            for model in [self.model_f, self.model_s]:
                model.initialize()
        else:  # restart
            self.model_f = self.restart_data['models'][0]
            self.model_s = self.restart_data['models'][1]
        for model in [self.model_f, self.model_s]:
            self.components += [model]

    def lop_f(self, dx):
        self.dxtemp.set_interface_data(dx.flatten())
        return self.model_f.predict(self.dxtemp).get_interface_data()

    def lop_s(self, dy):
        self.dytemp.set_interface_data(dy.flatten())
        return self.model_s.predict(self.dytemp).get_interface_data()

    def identity_matvec(self, v):
        return v

    def callback(self, rk):
        pass

    def solve_solution_step(self):
        iu = LinearOperator((self.u, self.u), self.identity_matvec)
        iw = LinearOperator((self.w, self.w), self.identity_matvec)
        dx = self.x.copy()
        dy = self.y.copy()
        # initial value
        self.x = self.predictor.predict(self.x)
        # first coupling iteration
        self.y = yt = self.solver_wrappers[0].solve_solution_step(self.x.copy()).copy()
        self.model_f.add(self.x.copy(), yt.copy())
        xt = self.solver_wrappers[1].solve_solution_step(self.y.copy()).copy()
        r = xt - self.x
        self.model_s.add(self.y.copy(), xt)
        self.finalize_iteration(r)
        # coupling iteration loop
        while not self.convergence_criterion.is_satisfied():
            if not self.model_s.is_ready() or not self.model_f.is_ready:
                dx = self.omega * r
            else:
                mf = LinearOperator((self.w, self.u), self.lop_f)
                ms = LinearOperator((self.u, self.w), self.lop_s)
                a = iu - ms @ mf
                b = (xt - self.x).get_interface_data() + ms @ (yt - self.y).get_interface_data()
                dx_sol, exitcode = gmres(a, b, tol=self.rtol, atol=self.atol, maxiter=20, callback=self.callback)
                if exitcode != 0:
                    RuntimeError('GMRES failed')
                dx.set_interface_data(dx_sol)
            self.x += dx
            yt = self.solver_wrappers[0].solve_solution_step(self.x.copy()).copy()
            self.model_f.add(self.x.copy(), yt.copy())
            if not self.model_s.is_ready() or not self.model_f.is_ready:
                dy = yt - self.y
            else:
                a = iw - mf @ ms
                b = (yt - self.y).get_interface_data() + mf @ (xt - self.x).get_interface_data()
                dy_sol, exitcode = gmres(a, b, tol=self.rtol, atol=self.atol, maxiter=20, callback=self.callback)
                if exitcode != 0:
                    RuntimeError('GMRES failed')
                dy.set_interface_data(dy_sol)
            self.y += dy
            xt = self.solver_wrappers[1].solve_solution_step(self.y.copy()).copy()
            r = xt - self.x
            self.model_s.add(self.y.copy(), xt)
            self.finalize_iteration(r)

    def check_restart_data(self, restart_data):
        for model in ['model_f', 'model_s']:
            model_original = self.parameters['settings'][model]['type']
            model_new = restart_data['parameters']['settings'][model]['type']
            if model_original != model_new:
                raise ValueError(f'Restart not possible because {model} type changed:'
                                 f'\n\toriginal: {model_original}\n\tnew: {model_new}')

    def add_restart_data(self, restart_data):
        return restart_data.update({'models': [self.model_f, self.model_s]})
