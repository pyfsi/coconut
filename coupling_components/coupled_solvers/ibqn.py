from coconut.coupling_components.tools import CreateInstance
from coconut.coupling_components.coupled_solvers.gauss_seidel import CoupledSolverGaussSeidel

from scipy.sparse.linalg import gmres, LinearOperator


def Create(parameters):
    return CoupledSolverIBQN(parameters)


class CoupledSolverIBQN(CoupledSolverGaussSeidel):
    def __init__(self, parameters):
        super().__init__(parameters)

        self.model_f = CreateInstance(self.parameters["settings"]["model_f"])
        self.model_s = CreateInstance(self.parameters["settings"]["model_s"])
        self.omega = self.settings["omega"].GetDouble()
        self.atol = self.settings["absolute_tolerance_gmres"].GetDouble()
        self.rtol = self.settings["relative_tolerance_gmres"].GetDouble()

        self.xtemp = self.ytemp = None
        self.dxtemp = self.dytemp = None
        self.u = self.w = None
        self.ready = None

    def initialize(self):
        super().initialize()

        self.dxtemp = self.x.deepcopy()
        self.dytemp = self.y.deepcopy()
        self.u = self.x.GetNumpyArray().shape[0]
        self.w = self.y.GetNumpyArray().shape[0]
        self.ready = False
        self.model_f.size_in = self.model_s.size_out = self.u
        self.model_f.size_out = self.model_s.size_in = self.w
        self.model_f.out = self.y.deepcopy()
        self.model_s.out = self.x.deepcopy()
        models = [self.model_f, self.model_s]
        for model in models:
            model.initialize()
            self.components += [model]

    def lop_f(self, dx):
        self.dxtemp.SetNumpyArray(dx.flatten())
        return self.model_f.Predict(self.dxtemp).GetNumpyArray()

    def lop_s(self, dy):
        self.dytemp.SetNumpyArray(dy.flatten())
        return self.model_s.Predict(self.dytemp).GetNumpyArray()

    def identity_matvec(self, v):
        return v

    def callback(self, rk):
        pass

    def solve_solution_step(self):
        iu = LinearOperator((self.u, self.u), self.identity_matvec)
        iw = LinearOperator((self.w, self.w), self.identity_matvec)
        dx = self.x.deepcopy()
        dy = self.y.deepcopy()
        # Initial value
        self.x = self.predictor.Predict(self.x)
        # First coupling iteration
        yt = self.solver_wrappers[0].solve_solution_step(self.x)
        self.model_f.Add(self.x, yt)
        self.y = yt
        xt = self.solver_wrappers[1].solve_solution_step(self.y)
        self.model_s.Add(self.y, xt)
        r = xt - self.x
        self.finalize_Iteration(r)
        # Coupling iteration loop
        while not self.convergence_criterion.is_satisfied():
            if not self.model_s.IsReady() or not self.model_f.IsReady:
                dx = self.omega * r
            else:
                mf = LinearOperator((self.w, self.u), self.lop_f)
                ms = LinearOperator((self.u, self.w), self.lop_s)
                a = iu - ms @ mf
                b = (xt - self.x).GetNumpyArray() + ms @ (yt - self.y).GetNumpyArray()
                dx_sol, exitcode = gmres(a, b, tol=self.rtol, atol=self.atol, maxiter=20, callback=self.callback)
                if exitcode != 0:
                    RuntimeError("GMRES failed")
                dx.SetNumpyArray(dx_sol)
            self.x += dx
            yt = self.solver_wrappers[0].solve_solution_step(self.x)
            self.model_f.Add(self.x, yt)
            if not self.model_s.IsReady() or not self.model_f.IsReady:
                dy = yt - self.y
            else:
                a = iw - mf @ ms
                b = (yt - self.y).GetNumpyArray() + mf @ (xt - self.x).GetNumpyArray()
                dy_sol, exitcode = gmres(a, b, tol=self.rtol, atol=self.atol, maxiter=20, callback=self.callback)
                if exitcode != 0:
                    RuntimeError("GMRES failed")
                dy.SetNumpyArray(dy_sol)
            self.y += dy
            xt = self.solver_wrappers[1].solve_solution_step(self.y)
            self.model_s.Add(self.y, xt)
            r = xt - self.x
            self.finalize_Iteration(r)
