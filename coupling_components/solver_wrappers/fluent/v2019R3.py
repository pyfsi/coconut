from coconut.coupling_components.solver_wrappers.fluent.fluent import SolverWrapperFluent
from coconut import tools


def create(parameters):
    return SolverWrapperFluent2019R3(parameters)


class SolverWrapperFluent2019R3(SolverWrapperFluent):

    def __init__(self, parameters):
        super().__init__(parameters)
        self.env = tools.get_solver_env(__name__, self.dir_cfd)
        self.check_software()

    def set_fluent_version(self):
        self.version = '2019R3'
        self.version_bis = '19.5.0'
