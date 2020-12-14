from coconut.coupling_components.convergence_criteria.combined import ConvergenceCriterionCombined


def create(parameters):
    return ConvergenceCriterionAnd(parameters)


class ConvergenceCriterionAnd(ConvergenceCriterionCombined):
    def IsSatisfied(self): #TODO: change to lower case
        is_satisfied = True
        for convergence_criterion in self.convergence_criteria:
            is_satisfied = is_satisfied and convergence_criterion.IsSatisfied()

        return is_satisfied
