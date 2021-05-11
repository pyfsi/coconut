from coconut.coupling_components.solver_wrappers.fluent.fluent import SolverWrapperFluent


def create(parameters):
    return SolverWrapperFluent2019R1(parameters)


class SolverWrapperFluent2019R1(SolverWrapperFluent):
    version = '2019R1'
    version_bis = '19.3.0'

    def __init__(self, parameters):
        super().__init__(parameters)
