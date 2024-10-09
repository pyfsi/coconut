from coconut.coupling_components.solver_wrappers.pc_fluent.fluent import SolverWrapperPCFluent
from coconut import tools


def create(parameters):
    return SolverWrapperPCFluent2023R1(parameters)


class SolverWrapperPCFluent2023R1(SolverWrapperPCFluent):
    version = '2023R1'
    version_bis = '23.1.0'

    def __init__(self, parameters):
        super().__init__(parameters)
        self.env = tools.get_solver_env(__name__, self.dir_cfd)
        self.check_software()
