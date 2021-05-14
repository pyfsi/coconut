from coconut.coupling_components.component import Component


def create(parameters):
    return ConvergenceCriterionIterationLimit(parameters)


class ConvergenceCriterionIterationLimit(Component):
    def __init__(self, parameters):
        super().__init__()

        settings = parameters['settings']
        self.maximum = settings['maximum']

        self.iteration = 0

    def initialize_solution_step(self):
        super().initialize_solution_step()

        self.iteration = 0

    def update(self, _unused):
        self.iteration += 1

    def is_satisfied(self):
        return self.iteration >= self.maximum
