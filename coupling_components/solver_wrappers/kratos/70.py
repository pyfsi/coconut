from coconut.coupling_components.component import Component


def Create(parameters):
    return SolverWrapperKratos70(parameters)


class SolverWrapperKratos70(Component):
    def __init__(self, parameters):
        super().__init__()

        self.settings = parameters["settings"]
