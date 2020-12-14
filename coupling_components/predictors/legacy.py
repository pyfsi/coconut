from coconut.coupling_components.predictors.predictor import Predictor


def create(parameters):
    return PredictorLegacy(parameters)


# Class PredictorQuadratic: Quadratic extrapolation based on the last three time steps, assuming constant time step size.
class PredictorLegacy(Predictor):
    def __init__(self, _unused):
        super().__init__(_unused)

        self.order = 2

    def predict(self, x):
        if len(self.dataprev) < 3:
            return self.linear(x)
        else:
            return self.legacy(x)
