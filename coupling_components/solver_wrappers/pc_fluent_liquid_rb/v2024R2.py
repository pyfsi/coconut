from coconut.coupling_components.solver_wrappers.pc_fluent_liquid_rb.fluent import SolverWrapperPCFluentLiquidRB
from coconut import tools


def create(parameters):
    return SolverWrapperPCFluent2024R2LiquidRB(parameters)


class SolverWrapperPCFluent2024R2LiquidRB(SolverWrapperPCFluentLiquidRB):
    version = '2024R2'
    version_bis = '24.2.0'

    def __init__(self, parameters):
        super().__init__(parameters)
        self.env = tools.get_solver_env(__name__, self.dir_cfd)
        self.check_software()
