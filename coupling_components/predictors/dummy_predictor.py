from coconut.coupling_components.component import Component


def create(parameters):
    return PredictorDummy(parameters)


class PredictorDummy(Component):
    def __init__(self, _):
        super().__init__()

    def predict(self, x):
        raise NotImplementedError('predict(x) called for PredictorDummy, use another predictor')

    def update(self, x):
        pass
