from coconut.coupling_components.component import CoSimulationComponent


def Create(parameters):
    return SolverWrapperKratos70(parameters)


class SolverWrapperKratos70(CoSimulationComponent):
    def __init__(self, parameters):
        super().__init__()

        self.settings = parameters["settings"]
