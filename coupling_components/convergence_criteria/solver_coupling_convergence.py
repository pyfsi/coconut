from coconut.coupling_components.component import Component
from coconut.coupling_components.solver_wrappers.solver_wrapper import SolverWrapper


def create(parameters):
    return ConvergenceCriterionSolverCouplingConvergence(parameters)


class ConvergenceCriterionSolverCouplingConvergence(Component):
    def __init__(self, parameters):
        super().__init__()

        settings = parameters['settings']
        self.solver_wrapper_index = settings['solver_index']
        self.solver_wrapper = None

    def initialize(self, solver_wrappers):
        super().initialize()

        if not (isinstance(solver_wrappers, list) and all(isinstance(s, SolverWrapper) for s in solver_wrappers)):
            raise ValueError(f'{self.__class__.__name__} has to be initialized with a list of solver wrappers')
        self.solver_wrapper = solver_wrappers[self.solver_wrapper_index]

        if not self.solver_wrapper.check_coupling_convergence_possible:
            raise ValueError(f'{self.solver_wrapper.__class__.__name__} does not allow to check coupling convergence')

        self.solver_wrapper.check_coupling_convergence = True

    def update(self, r):
        pass

    def is_satisfied(self):
        return self.solver_wrapper.coupling_convergence
