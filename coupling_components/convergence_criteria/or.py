from coconut.coupling_components.convergence_criteria.combined import ConvergenceCriterionCombined


def Create(parameters):
    return ConvergenceCriterionOr(parameters)


class ConvergenceCriterionOr(ConvergenceCriterionCombined):
    def IsSatisfied(self):
        is_satisfied = False
        for convergence_criterion in self.convergence_criteria:
            is_satisfied = is_satisfied or convergence_criterion.IsSatisfied()

        return is_satisfied
