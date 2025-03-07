from coconut.coupling_components.solver_wrappers.abaqus_cse.abaqus_cse import SolverWrapperAbaqusCSE
from coconut import tools


def create(parameters):
    return SolverWrapperAbaqusCSE2024(parameters)


class SolverWrapperAbaqusCSE2024(SolverWrapperAbaqusCSE):
    version = '2024'

    def __init__(self, parameters):
        super().__init__(parameters)
        self.env = tools.get_solver_env(__name__, self.dir_csm)
        self.check_software()
