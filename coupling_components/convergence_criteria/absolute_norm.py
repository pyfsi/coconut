from coconut.coupling_components.component import Component

import numpy as np


def Create(parameters):
    return ConvergenceCriterionAbsoluteNorm(parameters)


class ConvergenceCriterionAbsoluteNorm(Component):
    def __init__(self, parameters):
        super().__init__()

        settings = parameters["settings"]
        self.tolerance = settings["tolerance"].GetDouble()
        self.order = settings["order"].GetInt()

        self.last_norm = 0.0
        self.is_updated = False

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()

        self.last_norm = 0.0
        self.is_updated = False

    def Update(self, r):
        self.last_norm = np.linalg.norm(r.GetNumpyArray(), self.order)
        self.is_updated = True

    def IsSatisfied(self):
        # tools.PrintInfo("Norm: " + str(self.last_norm))
        if not self.is_updated:
            return False
        else:
            return self.last_norm < self.tolerance
