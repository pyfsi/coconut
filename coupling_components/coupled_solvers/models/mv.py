from coconut.coupling_components.component import Component
from coconut import tools

import numpy as np


def create(parameters):
    return ModelMV(parameters)


class ModelMV(Component):
    def __init__(self, parameters):
        super().__init__()

        self.settings = parameters["settings"]
        self.min_significant = self.settings["min_significant"]

        self.size_in = None
        self.size_out = None
        self.out = None  # interface of output
        self.added = False
        self.rref = None
        self.xtref = None
        self.v = None
        self.w = None
        self.ncurr = None
        self.nprev = None

    def initialize(self):
        super().initialize()

        self.v = np.empty((self.size_in, 0))
        self.w = np.empty((self.size_out, 0))
        self.nprev = np.zeros((self.size_out, self.size_in))

    def filter(self):
        if self.v.shape[1] == 0:
            raise RuntimeError("No information to filter")
        # remove columns resulting in small diagonal elements in R
        singular = True
        while singular and self.v.shape[1]:
            rr = np.linalg.qr(self.v, mode='r')
            diag = np.diagonal(rr)
            m = min(abs(diag))
            if m < self.min_significant:
                i = np.argmin(abs(diag))
                tools.print_info("Removing column " + str(i) + ": " + str(m) + " < minsignificant", layout='warning')
                self.v = np.delete(self.v, i, 1)
                self.w = np.delete(self.w, i, 1)
            else:
                singular = False
        # remove columns if number of columns exceeds number of rows
        if self.v.shape[0] < self.v.shape[1]:
            self.v = np.delete(self.v, -1, 1)
            self.w = np.delete(self.w, -1, 1)

    def predict(self, dr_in):
        dr = dr_in.get_interface_data().reshape(-1, 1)
        if self.ncurr is None:
            raise RuntimeError("No information to predict")
        # approximation for the inverse of the Jacobian from a multiple vector model
        dxt = self.ncurr @ dr
        dxt_out = self.out.copy()
        dxt_out.set_interface_data(dxt.flatten())
        return dxt_out

    def add(self, r_in, xt_in):
        r = r_in.get_interface_data().reshape(-1, 1)
        xt = xt_in.get_interface_data().reshape(-1, 1)
        if self.added:
            dr = r - self.rref
            dxt = xt - self.xtref
            # update V and W matrices
            self.v = np.hstack((dr, self.v))
            self.w = np.hstack((dxt, self.w))
            self.filter()
            # update of the matrix N
            self.ncurr = self.nprev + (self.w - self.nprev @ self.v) @ np.linalg.inv(self.v.T @ self.v) @ self.v.T
        else:
            self.added = True
        self.rref = r
        self.xtref = xt

    def is_ready(self):
        return self.ncurr is not None

    def initialize_solution_step(self):
        super().initialize_solution_step()

        self.rref = None
        self.xtref = None
        self.v = np.empty((self.size_in, 0))
        self.w = np.empty((self.size_out, 0))
        self.added = False

    def finalize_solution_step(self):
        super().finalize_solution_step()

        self.nprev = self.ncurr

    def filter_q(self, r_in):
        r = r_in.get_interface_data().reshape(-1, 1)
        r_out = r_in.copy()
        qt, *_ = np.linalg.qr(self.ncurr.T)
        q = qt[:, :np.linalg.matrix_rank(self.ncurr)]
        r = r - q @ (q.T @ r)
        r_out.set_interface_data(r.flatten())
        return r_out
