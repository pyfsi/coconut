from coconut.coupling_components.solver_wrappers.fluent.alm_fluent import SolverWrapperALMFluent
from coconut import tools


def create(parameters):
    return SolverWrapperALMFluent2023R1(parameters)


class SolverWrapperALMFluent2023R1(SolverWrapperALMFluent):
    version = '2023R1'
    version_bis = '23.1.0'

    def __init__(self, parameters):
        super().__init__(parameters)
        self.env = tools.get_solver_env(__name__, self.dir_cfd)
        self.check_software()
