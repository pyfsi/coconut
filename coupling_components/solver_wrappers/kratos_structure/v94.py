from coconut.coupling_components.solver_wrappers.kratos_structure.kratos_structure import SolverWrapperKratosStructure
from coconut import tools


def create(parameters):
    return SolverWrapperKratosStructure94(parameters)


class SolverWrapperKratosStructure94(SolverWrapperKratosStructure):
    version = '94'

    def __init__(self, parameters):
        super().__init__(parameters)
        self.env = tools.get_solver_env(__name__, self.working_directory)
        self.check_software()
