from coconut.coupling_components.component import Component
from coconut.coupling_components.tools import create_instance, PrintStructureInfo


class ConvergenceCriterionCombined(Component):
    def __init__(self, parameters):
        super().__init__()

        settings = parameters["settings"]
        self.convergence_criteria = []
        for criterion in settings["criteria_list"]:
            self.convergence_criteria.append(create_instance(criterion))


    def Initialize(self): #TODO: change to lower case
        super().Initialize()

        for convergence_criterion in self.convergence_criteria:
            convergence_criterion.Initialize()

    def Finalize(self): #TODO: change to lower case
        super().Finalize()

        for convergence_criterion in self.convergence_criteria:
            convergence_criterion.Finalize()

    def InitializeSolutionStep(self): #TODO: change to lower case
        super().InitializeSolutionStep()

        for convergence_criterion in self.convergence_criteria:
            convergence_criterion.InitializeSolutionStep()

    def FinalizeSolutionStep(self): #TODO: change to lower case
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
