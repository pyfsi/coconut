from coconut.coupling_components.solver_wrappers.fluent.fluent import SolverWrapperFluent


def create(parameters):
    return SolverWrapperFluent2020R1(parameters)


class SolverWrapperFluent2020R1(SolverWrapperFluent):
    version = '2020R1'
    version_bis = '20.1.0'

    def __init__(self, parameters):
        super().__init__(parameters)
