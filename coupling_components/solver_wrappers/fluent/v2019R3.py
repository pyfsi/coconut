from coconut.coupling_components.solver_wrappers.fluent.v2019R2 import SolverWrapperFluent2019R2


def create(parameters):
    return SolverWrapperFluent2019R3(parameters)


class SolverWrapperFluent2019R3(SolverWrapperFluent2019R2):

    def set_fluent_version(self):
        self.version = '2019R3'
        self.version_bis = '19.5.0'
