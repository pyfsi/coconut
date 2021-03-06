from coconut.coupling_components.solver_wrappers.fluent.v2019R1 import SolverWrapperFluent2019R1


def create(parameters):
    return SolverWrapperFluent2019R2(parameters)


class SolverWrapperFluent2019R2(SolverWrapperFluent2019R1):

    def set_fluent_version(self):
        self.version = '2019R2'
        self.version_bis = '19.4.0'
