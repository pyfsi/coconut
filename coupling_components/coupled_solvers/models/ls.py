from coconut.coupling_components.component import Component
from coconut.coupling_components import tools

import numpy as np
from scipy.linalg import solve_triangular


def Create(parameters):
    return ModelLS(parameters)


class ModelLS(Component):
    def __init__(self, parameters):
        super().__init__()

        self.settings = parameters["settings"]
        self.min_significant = self.settings["min_significant"].GetDouble()
        self.q = self.settings["q"].GetDouble()

        self.size_in = None
        self.size_out = None
        self.out = None  # Interface of output
        self.added = False
        self.rref = None
        self.xtref = None
        self.vcurr = None
        self.wcurr = None
        self.vprev = None
        self.wprev = None

    def Initialize(self):
        super().Initialize()

        self.vcurr = np.empty((self.size_in, 0))
        self.wcurr = np.empty((self.size_out, 0))
        self.vprev = [np.empty((self.size_in, 0))]
        self.wprev = [np.empty((self.size_out, 0))]

    def Filter(self):
        v = np.hstack((self.vcurr, np.hstack(self.vprev)))
        if not v.shape[1]:
            raise RuntimeError("No information to filter")
        # Remove columns resulting in small diagonal elements in R
        singular = True
        while singular and v.shape[1]:
            rr = np.linalg.qr(v, mode='r')
            diag = np.diagonal(rr)
            m = min(abs(diag))
            if m < self.min_significant:
                i = np.argmin(abs(diag))
                tools.Print("Removing column " + str(i) + ": " + str(m) + " < minsignificant", layout='warning')
                if i < self.vcurr.shape[1]:
                    self.vcurr = np.delete(self.vcurr, i, 1)
                    self.wcurr = np.delete(self.wcurr, i, 1)
                else:
                    num_columns = self.vcurr.shape[1]
                    j = -1
                    while i >= num_columns:
                        j += 1
                        num_columns += self.vprev[j].shape[1]
                    num_columns -= self.vprev[j].shape[1]
                    self.vprev[j] = np.delete(self.vprev[j], i - num_columns, 1)
                    self.wprev[j] = np.delete(self.wprev[j], i - num_columns, 1)
                v = np.hstack((self.vcurr, np.hstack(self.vprev)))
            else:
                singular = False
        # Remove columns if number of columns exceeds number of rows
        while v.shape[0] < v.shape[1]:
            if self.vcurr.shape[0] < self.vcurr.shape[1]:
                self.vcurr = np.delete(self.vcurr, -1, 1)
                self.wcurr = np.delete(self.wcurr, -1, 1)
            else:
                i = -1
                while self.vprev[i].shape[1] == 0:
                    i -= 1
                self.vprev[i] = np.delete(self.vprev[i], -1, 1)
                self.wprev[i] = np.delete(self.wprev[i], -1, 1)
            v = np.hstack((self.vcurr, np.hstack(self.vprev)))

    def Predict(self, dr_in):
        dr = dr_in.GetNumpyArray().reshape(-1, 1)
        v = np.hstack((self.vcurr, np.hstack(self.vprev)))
        w = np.hstack((self.wcurr, np.hstack(self.wprev)))
        if not v.shape[1]:
            raise RuntimeError("No information to predict")
        # Approximation for the inverse of the Jacobian from a least-squares model
        qq, rr = np.linalg.qr(v, mode='reduced')
        b = qq.T @ dr
        c = solve_triangular(rr, b)
        dxt = w @ c
        dxt_out = self.out.deepcopy()
        dxt_out.SetNumpyArray(dxt.flatten())
        return dxt_out

    def Add(self, r_in, xt_in):
        r = r_in.GetNumpyArray().reshape(-1, 1)
        xt = xt_in.GetNumpyArray().reshape(-1, 1)
        if self.added:
            dr = r - self.rref
            dxt = xt - self.xtref
            # Update V and W matrices
            self.vcurr = np.hstack((dr, self.vcurr))
            self.wcurr = np.hstack((dxt, self.wcurr))
            self.Filter()
        else:
            self.added = True
        self.rref = r
        self.xtref = xt

    def IsReady(self):
        v = np.hstack((self.vcurr, np.hstack(self.vprev)))
        return v.shape[1]

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()

        self.rref = None
        self.xtref = None
        self.vcurr = np.empty((self.size_in, 0))
        self.wcurr = np.empty((self.size_out, 0))
        self.added = False

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()

        if self.q > 0:
            self.vprev = [self.vcurr] + self.vprev
            self.wprev = [self.wcurr] + self.wprev
            if len(self.vprev) > self.q:
                self.vprev.pop()
                self.wprev.pop()

    def FilterQ(self, r_in):
        r = r_in.GetNumpyArray().reshape(-1, 1)
        r_out = r_in.deepcopy()
        v = np.hstack((self.vcurr, np.hstack(self.vprev)))
        if v.shape[1]:
            qq, *_ = np.linalg.qr(v, mode='reduced')
            r = r - qq @ (qq.T @ r)
            r_out.SetNumpyArray(r.flatten())
        return r_out