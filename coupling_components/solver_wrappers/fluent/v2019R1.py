from coconut.coupling_components.solver_wrappers.fluent.fluent import SolverWrapperFluent
from coconut import tools
import warnings


def create(parameters):
    return SolverWrapperFluent2019R1(parameters)


class SolverWrapperFluent2019R1(SolverWrapperFluent):

    def __init__(self, parameters):
        super().__init__(parameters)
        with warnings.catch_warnings():
            warnings.filterwarnings('always', category=DeprecationWarning)
            warnings.warn('SolverWrapperFluent2019R1 will no longer be maintained and tested', category=DeprecationWarning)
        self.env = tools.get_solver_env(__name__, self.dir_cfd)
        self.check_software()

    def set_fluent_version(self):
        self.version = '2019R1'
        self.version_bis = '19.3.0'
