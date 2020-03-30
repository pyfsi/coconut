from coconut.coupling_components.component import Component
from coconut.coupling_components import tools


def Create(parameters):
    return ConvergenceCriterionIterationLimit(parameters)


class ConvergenceCriterionIterationLimit(Component):
    def __init__(self, parameters):
        super().__init__()

        settings = parameters["settings"]
        self.maximum = settings["maximum"].GetInt()

        self.iteration = 0

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()

        self.iteration = 0

    def Update(self, _unused):
        self.iteration += 1

    def IsSatisfied(self):
        # tools.PrintInfo("Iteration: " + str(self.iteration))
        return self.iteration >= self.maximum
