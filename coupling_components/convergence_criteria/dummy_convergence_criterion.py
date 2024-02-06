from coconut.coupling_components.component import Component


def create(parameters):
    return ConvergenceCriterionDummy(parameters)


class ConvergenceCriterionDummy(Component):
    dummy = True

    def __init__(self, _):
        super().__init__()

    def update(self, r):
        pass

    def is_satisfied(self):
        raise NotImplementedError('is_satisfied() called for ConvergenceCriterionDummy,'
                                  ' use another convergence_criterion')
