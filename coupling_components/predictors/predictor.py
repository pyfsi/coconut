from coconut.coupling_components.component import Component


def create(parameters):
    return Predictor(parameters)


# base class for extrapolation based on the last two (linear), three (quadratic) or four (cubic) time steps,
# assuming constant time step size
class Predictor(Component):
    def __init__(self, _unused):
        super().__init__()

        self.updated = False
        self.dataprev = None
        self.order = None

    def initialize(self, x):
        super().initialize()

        self.dataprev = [x.get_interface_data()]

    def initialize_solution_step(self):
        super().initialize_solution_step()

        self.updated = False

    def finalize_solution_step(self):
        super().finalize_solution_step()
        if not self.updated:
            raise Exception('Not updated')

    def constant(self, x_in):
        x = x_in.copy()
        if not self.updated:
            y = self.dataprev[0]
            x.set_interface_data(y)
            return x
        else:
            raise Exception('Already updated')

    def linear(self, x_in):
        x = x_in.copy()
        if not self.updated:
            if len(self.dataprev) == 1:
                y = self.dataprev[0]
            else:
                y = 2 * self.dataprev[0] - self.dataprev[1]
            x.set_interface_data(y)
            return x
        else:
            raise Exception('Already updated')

    def quadratic(self, x_in):
        x = x_in.copy()
        if not self.updated:
            if len(self.dataprev) < 3:
                raise Exception('Not sufficient information for quadratic extrapolation')
            y = 3.0 * self.dataprev[0] - 3.0 * self.dataprev[1] + 1.0 * self.dataprev[2]
            x.set_interface_data(y)
            return x
        else:
            raise Exception('Already updated')

    def legacy(self, x_in):
        x = x_in.copy()
        if not self.updated:
            if len(self.dataprev) < 3:
                raise Exception('Not sufficient information for quadratic extrapolation')
            y = 2.5 * self.dataprev[0] - 2.0 * self.dataprev[1] + 0.5 * self.dataprev[2]
            x.set_interface_data(y)
            return x
        else:
            raise Exception('Already updated')

    def cubic(self, x_in):
        x = x_in.copy()
        if not self.updated:
            if len(self.dataprev) < 4:
                raise Exception('Not sufficient information for cubic extrapolation')
            y = 4.0 * self.dataprev[0] - 6.0 * self.dataprev[1] + 4.0 * self.dataprev[2] - 1.0 * self.dataprev[3]
            x.set_interface_data(y)
            return x
        else:
            raise Exception('Already updated')

    def predict(self, x):
        pass

    def update(self, x):
        if not self.updated:
            self.dataprev = [x.get_interface_data()] + self.dataprev
            if len(self.dataprev) > self.order + 1:
                self.dataprev.pop()
            self.updated = True
        else:
            raise Exception('Already updated')

    def restart(self, restart_data):
        self.dataprev = restart_data['dataprev']
        while len(self.dataprev) > self.order + 1:
            self.dataprev.pop()

    def check_restart_data(self, restart_data):
        pass

    def save_restart_data(self):
        return {'dataprev': self.dataprev}
