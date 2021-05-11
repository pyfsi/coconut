from coconut.coupling_components.solver_wrappers.fluent.fluent import SolverWrapperFluent


def create(parameters):
    return SolverWrapperFluent2019R2(parameters)


class SolverWrapperFluent2019R2(SolverWrapperFluent):
    version = '2019R2'
    version_bis = '19.4.0'

    def __init__(self, parameters):
        super().__init__(parameters)
