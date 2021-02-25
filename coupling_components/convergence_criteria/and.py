from coconut.coupling_components.convergence_criteria.combined import ConvergenceCriterionCombined


def create(parameters):
    return ConvergenceCriterionAnd(parameters)


class ConvergenceCriterionAnd(ConvergenceCriterionCombined):
    def is_satisfied(self):
        is_satisfied = True
        for convergence_criterion in self.convergence_criteria:
            is_satisfied = is_satisfied and convergence_criterion.is_satisfied()

        return is_satisfied
