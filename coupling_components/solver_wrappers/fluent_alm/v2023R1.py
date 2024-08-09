from coconut.coupling_components.solver_wrappers.fluent_alm.fluent_alm import SolverWrapperFluentALM
from coconut import tools


def create(parameters):
    return SolverWrapperFluentALM2023R1(parameters)


class SolverWrapperFluentALM2023R1(SolverWrapperFluentALM):
    version = '2023R1'
    version_bis = '23.1.0'

    def __init__(self, parameters):
        super().__init__(parameters)
        self.env = tools.get_solver_env(__name__, self.dir_cfd)
        self.check_software()
