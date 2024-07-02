from coconut.coupling_components.solver_wrappers.abaqus.abaqus import SolverWrapperAbaqus
from coconut import tools


def create(parameters):
    return SolverWrapperAbaqus2024(parameters)


class SolverWrapperAbaqus2024(SolverWrapperAbaqus):
    version = '2024'

    def __init__(self, parameters):
        super().__init__(parameters)
        self.env = tools.get_solver_env(__name__, self.dir_csm)
        self.check_software()

        # environment file parameters
        self.link_sl = '-Wl,--add-needed'
        self.link_exe = '-Wl,--add-needed'
