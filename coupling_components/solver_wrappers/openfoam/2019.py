from coconut.coupling_components.component import Component


def create(parameters):
    return SolverWrapperOpenFOAM2019(parameters)


class SolverWrapperOpenFOAM2019(Component):
    def __init__(self, parameters):
        super().__init__()

        self.settings = parameters["settings"]
