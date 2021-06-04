from coconut.coupling_components.component import Component


def create(parameters):
    return ConvergenceCriterionAbsoluteNorm(parameters)


class ConvergenceCriterionAbsoluteNorm(Component):
    def __init__(self, parameters):
        super().__init__()

        settings = parameters['settings']
        self.tolerance = settings['tolerance']
        self.order = settings['order']

        self.last_norm = 0
        self.is_updated = False

    def initialize_solution_step(self):
        super().initialize_solution_step()

        self.last_norm = 0
        self.is_updated = False

    def update(self, r):
        self.last_norm = r.norm(order=self.order)
        self.is_updated = True

    def is_satisfied(self):
        if not self.is_updated:
            return False
        else:
            return self.last_norm < self.tolerance
