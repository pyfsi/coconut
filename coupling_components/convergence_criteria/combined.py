from coconut.coupling_components.component import Component
from coconut.coupling_components.tools import CreateInstance, PrintStructureInfo


class ConvergenceCriterionCombined(Component):
    def __init__(self, parameters):
        super().__init__()

        settings = parameters["settings"]
        self.convergence_criteria = []
        index = 0
        while True:
            key = "convergence_criterion" + str(index)
            if key in settings.keys():
                self.convergence_criteria.append(CreateInstance(settings[key]))
                index += 1
            else:
                break

    def Initialize(self):
        super().Initialize()

        for convergence_criterion in self.convergence_criteria:
            convergence_criterion.Initialize()

    def Finalize(self):
        super().Finalize()

        for convergence_criterion in self.convergence_criteria:
            convergence_criterion.Finalize()

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()

        for convergence_criterion in self.convergence_criteria:
            convergence_criterion.InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()

        for convergence_criterion in self.convergence_criteria:
            convergence_criterion.FinalizeSolutionStep()

    def OutputSolutionStep(self):
        super().OutputSolutionStep()

        for convergence_criterion in self.convergence_criteria:
            convergence_criterion.OutputSolutionStep()

    def Check(self):
        super().Check()

        for convergence_criterion in self.convergence_criteria:
            convergence_criterion.Check()

    def PrintInfo(self, pre):
        super().PrintInfo(pre)

        PrintStructureInfo(pre, self.convergence_criteria)

    def Update(self, r):
        for convergence_criterion in self.convergence_criteria:
            convergence_criterion.Update(r)
