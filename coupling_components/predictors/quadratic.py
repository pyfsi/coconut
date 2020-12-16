from coconut.coupling_components.predictors.predictor import Predictor


def create(parameters):
    return PredictorQuadratic(parameters)


# quadratic extrapolation based on the last three time steps, assuming constant time step size
class PredictorQuadratic(Predictor):
    def __init__(self, _unused):
        super().__init__(_unused)

        self.order = 2

    def predict(self, x):
        if len(self.dataprev) < 3:
            return self.linear(x)
        else:
            return self.quadratic(x)
