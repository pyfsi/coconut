from coconut.coupling_components.predictors.predictor import Predictor


def Create(parameters):
    return PredictorConstant(parameters)


# Class PredictorLinear: Linear extrapolation based on the last two time steps, assuming constant time step size.
class PredictorConstant(Predictor):
    def __init__(self, _unused):
        super().__init__(_unused)

        self.order = 0

    def Predict(self, x):
        return self.Constant(x)
