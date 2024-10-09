from coconut.coupling_components.solver_wrappers.cht_fluent.fluent import SolverWrapperCHTFluent
from coconut import tools


def create(parameters):
    return SolverWrapperCHTFluent2023R1(parameters)


class SolverWrapperCHTFluent2023R1(SolverWrapperCHTFluent):
    version = '2023R1'
    version_bis = '23.1.0'

    def __init__(self, parameters):
        super().__init__(parameters)
        self.env = tools.get_solver_env(__name__, self.dir_cfd)
        self.check_software()
