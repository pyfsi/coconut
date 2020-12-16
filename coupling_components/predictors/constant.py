from coconut.coupling_components.predictors.predictor import Predictor


def create(parameters):
    return PredictorConstant(parameters)


# constant extrapolation based on the last time step
class PredictorConstant(Predictor):
    def __init__(self, _unused):
        super().__init__(_unused)

        self.order = 0

    def predict(self, x):
        return self.constant(x)
