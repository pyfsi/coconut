from coconut.coupling_components.predictors.predictor import Predictor


def Create(parameters):
    return PredictorConstant(parameters)


# Class PredictorConstant: Constant extrapolation based on the last time step.
class PredictorConstant(Predictor):
    def __init__(self, _unused):
        super().__init__(_unused)

        self.order = 0

    def Predict(self, x):
        return self.Constant(x)
