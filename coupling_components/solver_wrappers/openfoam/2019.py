from coconut.coupling_components.component import CoSimulationComponent


def Create(parameters):
    return SolverWrapperOpenFOAM2019(parameters)


class SolverWrapperOpenFOAM2019(CoSimulationComponent):
    def __init__(self, parameters):
        super().__init__()

        self.settings = parameters["settings"]
