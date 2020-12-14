from coconut.coupling_components.convergence_criteria.combined import ConvergenceCriterionCombined


def create(parameters):
    return ConvergenceCriterionOr(parameters)


class ConvergenceCriterionOr(ConvergenceCriterionCombined):
    def is_satisfied(self):
        is_satisfied = False
        for convergence_criterion in self.convergence_criteria:
            is_satisfied = is_satisfied or convergence_criterion.is_satisfied()

        return is_satisfied
