from coconut.coupling_components.solver_wrappers.openfoam_extend.openfoam_extend import SolverWrapperOpenFOAMExtend
from coconut import tools
from coconut.coupling_components.solver_wrappers.openfoam import openfoam_io as of_io

from subprocess import check_call
from os.path import join


def create(parameters):
    return SolverWrapperOpenFOAMExtend41(parameters)


class SolverWrapperOpenFOAMExtend41(SolverWrapperOpenFOAMExtend):
    version = '4.1'

    def __init__(self, parameters):
        super().__init__(parameters)
        self.env = tools.get_solver_env(__name__, self.working_directory)

        # compile adapted openfoam software
        self.compile_adapted_openfoam_extend_solver()

        # check that the correct software is available
        self.check_software()





