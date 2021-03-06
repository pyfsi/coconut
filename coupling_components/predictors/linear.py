from coconut.coupling_components.predictors.predictor import Predictor


def create(parameters):
    return PredictorLinear(parameters)


# linear extrapolation based on the last two time steps, assuming constant time step size
class PredictorLinear(Predictor):
    def __init__(self, _unused):
        super().__init__(_unused)

        self.order = 1

    def predict(self, x):
        return self.linear(x)
