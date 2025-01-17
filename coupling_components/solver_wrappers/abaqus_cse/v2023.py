from coconut.coupling_components.solver_wrappers.abaqus_cse.abaqus_cse import SolverWrapperAbaqusCSE
from coconut import tools


def create(parameters):
    return SolverWrapperAbaqusCSE2023(parameters)


class SolverWrapperAbaqusCSE2023(SolverWrapperAbaqusCSE):
    version = '2023'

    def __init__(self, parameters):
        super().__init__(parameters)
        self.env = tools.get_solver_env(__name__, self.dir_csm)
        self.check_software()
