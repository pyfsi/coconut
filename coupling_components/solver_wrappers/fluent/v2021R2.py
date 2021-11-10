from coconut.coupling_components.solver_wrappers.fluent.fluent import SolverWrapperFluent
from coconut import tools


def create(parameters):
    return SolverWrapperFluent2021R2(parameters)


class SolverWrapperFluent2021R2(SolverWrapperFluent):
    version = '2021R2'
    version_bis = '21.2.0'

    def __init__(self, parameters):
        super().__init__(parameters)
        self.env = tools.get_solver_env(__name__, self.dir_cfd)
        self.check_software()
