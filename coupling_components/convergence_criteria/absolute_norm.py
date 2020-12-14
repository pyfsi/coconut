from coconut.coupling_components.component import Component

import numpy as np


def create(parameters):
    return ConvergenceCriterionAbsoluteNorm(parameters)


class ConvergenceCriterionAbsoluteNorm(Component):
    def __init__(self, parameters):
        super().__init__()

        settings = parameters["settings"]
        self.tolerance = settings["tolerance"]
        self.order = settings["order"]

        self.last_norm = 0.0
        self.is_updated = False

    def InitializeSolutionStep(self): #TODO: change to lower case
        super().InitializeSolutionStep()

        self.last_norm = 0.0
        self.is_updated = False

    def Update(self, r): #TODO: change to lower case
        self.last_norm = np.linalg.norm(r.get_interface_data(), self.order)
        self.is_updated = True

    def IsSatisfied(self): #TODO: change to lower case
        # tools.PrintInfo("Norm: " + str(self.last_norm))
        if not self.is_updated:
            return False
        else:
            return self.last_norm < self.tolerance
