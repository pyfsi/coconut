from coconut.coupling_components.solver_wrappers.fluent.fluent import SolverWrapperFluent


def create(parameters):
    return SolverWrapperFluent2019R3(parameters)


class SolverWrapperFluent2019R3(SolverWrapperFluent):
    version = '2019R3'
    version_bis = '19.5.0'

    def __init__(self, parameters):
        super().__init__(parameters)
