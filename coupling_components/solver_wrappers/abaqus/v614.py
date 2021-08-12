from coconut.coupling_components.solver_wrappers.abaqus.abaqus import SolverWrapperAbaqus
from coconut import tools


def create(parameters):
    return SolverWrapperAbaqus614(parameters)


class SolverWrapperAbaqus614(SolverWrapperAbaqus):
    version = '6.14'

    def __init__(self, parameters):
        super().__init__(parameters)
        self.env = tools.get_solver_env(__name__, self.dir_csm)
        self.check_software()

        # environment file parameters
        self.link_sl = '-i-dynamic'
