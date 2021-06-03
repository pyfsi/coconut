from coconut.coupling_components.component import Component

import numpy as np


def create(parameters):
    return ConvergenceCriterionRelativeNorm(parameters)


class ConvergenceCriterionRelativeNorm(Component):
    def __init__(self, parameters):
        super().__init__()

        settings = parameters['settings']
        self.tolerance = settings['tolerance']
        self.order = settings['order']

        self.initial_norm = 0
        self.last_norm = 0
        self.is_initial_norm_set = False

    def initialize_solution_step(self):
        super().initialize_solution_step()

        self.initial_norm = 0
        self.last_norm = 0
        self.is_initial_norm_set = False

    def update(self, r):
        self.last_norm = r.norm(order=self.order)
        if not self.is_initial_norm_set:
            self.initial_norm = self.last_norm
            self.is_initial_norm_set = True
            if self.initial_norm < np.finfo(type(self.initial_norm)).eps:
                raise Exception('Initial norm is too small')

    def is_satisfied(self):
        if not self.is_initial_norm_set:
            return False
        else:
            return self.last_norm / self.initial_norm < self.tolerance
