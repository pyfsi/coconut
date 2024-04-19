from coconut.coupling_components.solver_wrappers.fluent.pc_fluent_iso import SolverWrapperFluent
from coconut import tools


def create(parameters):
    return SolverWrapperFluent2023R1(parameters)


class SolverWrapperFluent2023R1(SolverWrapperFluent):
    version = '2023R1'
    version_bis = '23.1.0'

    def __init__(self, parameters):
        super().__init__(parameters)
        self.env = tools.get_solver_env(__name__, self.dir_cfd)
        self.check_software()
