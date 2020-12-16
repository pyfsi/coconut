from coconut.coupling_components.component import Component
from coconut.coupling_components.tools import create_instance, print_components_info


class ConvergenceCriterionCombined(Component):
    def __init__(self, parameters):
        super().__init__()

        settings = parameters["settings"]
        self.convergence_criteria = []
        for criterion in settings["criteria_list"]:
            self.convergence_criteria.append(create_instance(criterion))

    def initialize(self):
        super().initialize()

        for convergence_criterion in self.convergence_criteria:
            convergence_criterion.initialize()

    def finalize(self):
        super().finalize()

        for convergence_criterion in self.convergence_criteria:
            convergence_criterion.finalize()

    def initialize_solution_step(self):
        super().initialize_solution_step()

        for convergence_criterion in self.convergence_criteria:
            convergence_criterion.initialize_solution_step()

    def finalize_solution_step(self):
        super().finalize_solution_step()

        for convergence_criterion in self.convergence_criteria:
            convergence_criterion.finalize_solution_step()

    def output_solution_step(self):
        super().output_solution_step()

        for convergence_criterion in self.convergence_criteria:
            convergence_criterion.output_solution_step()

    def check(self):
        super().check()

        for convergence_criterion in self.convergence_criteria:
            convergence_criterion.check()

    def print_components_info(self, pre):
        super().print_components_info(pre)

        print_components_info(pre, self.convergence_criteria)

    def update(self, r):
        for convergence_criterion in self.convergence_criteria:
            convergence_criterion.update(r)
