from coconut.coupling_components.solver_wrappers.abaqus.abaqus_line_load import SolverWrapperAbaqusLineLoad
from coconut import tools


def create(parameters):
    return SolverWrapperAbaqusLineLoad2024(parameters)


class SolverWrapperAbaqusLineLoad2024(SolverWrapperAbaqusLineLoad):
    version = '2024'

    def __init__(self, parameters):
        super().__init__(parameters)
        self.env = tools.get_solver_env(__name__, self.dir_csm)
        self.check_software()

        # environment file parameters
        self.link_sl = '-Wl,--add-needed'
        self.link_exe = '-Wl,--add-needed'
