from coconut.coupling_components.solver_wrappers.openfoam.openfoam import SolverWrapperOpenFOAM
from coconut import tools


def create(parameters):
    return SolverWrapperOpenFOAM41(parameters)


class SolverWrapperOpenFOAM41(SolverWrapperOpenFOAM):
    version = '4.1'

    def __init__(self, parameters):
        super().__init__(parameters)
        self.env = tools.get_solver_env(__name__, self.working_directory)

        # compile adapted openfoam software
        self.compile_adapted_openfoam_solver()

        # check that the correct software is available
        self.check_software()
