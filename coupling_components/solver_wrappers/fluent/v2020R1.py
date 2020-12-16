from coconut.coupling_components.solver_wrappers.fluent.v2019R3 import SolverWrapperFluent2019R3


def create(parameters):
    return SolverWrapperFluent2020R1(parameters)


class SolverWrapperFluent2020R1(SolverWrapperFluent2019R3):

    def set_fluent_version(self):
        self.version = '2020R1'
        self.version_bis = '20.1.0'
