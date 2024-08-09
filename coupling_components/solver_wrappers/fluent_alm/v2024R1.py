from coconut.coupling_components.solver_wrappers.fluent_alm.fluent_alm import SolverWrapperFluentALM
from coconut import tools


def create(parameters):
    return SolverWrapperFluentALM2024R1(parameters)


class SolverWrapperFluentALM2024R1(SolverWrapperFluentALM):
    version = '2024R1'
    version_bis = '24.1.0'

    def __init__(self, parameters):
        super().__init__(parameters)
        self.env = tools.get_solver_env(__name__, self.dir_cfd)
        self.check_software()
