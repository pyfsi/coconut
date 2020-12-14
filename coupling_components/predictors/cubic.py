from coconut.coupling_components.predictors.predictor import Predictor


def create(parameters):
    return PredictorCubic(parameters)


# Class PredictorCubic: Cubic extrapolation based on the last four time steps, assuming constant time step size.
class PredictorCubic(Predictor):
    def __init__(self, _unused):
        super().__init__(_unused)

        self.order = 3

    def predict(self, x):
        if len(self.dataprev) < 3:
            return self.linear(x)
        elif len(self.dataprev) == 3:
            return self.quadratic(x)
        else:
            return self.cubic(x)
