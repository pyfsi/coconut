from coconut.coupling_components.component import Component
from coconut import tools


def create(parameters):
    return Surrogate(parameters)


class Surrogate(Component):
    def __init__(self, parameters):
        super().__init__()

        self.parameters = parameters
        if 'settings' not in self.parameters:
            self.parameters['settings'] = {}
        self.settings = parameters.get('settings')

        # read parameters
        self.predict_change = self.settings.get('predict_change', True)  # use surrogate change rather than value

        self.updated = False
        self.surrogate_updated = False
        self.dataprev = None
        self.surrogate_data = None
        self.surrogate_dataprev = None
        self.order = None

    def initialize(self, x):
        super().initialize()

        self.dataprev = x.get_interface_data()
        self.surrogate_data = x.get_interface_data()

    def initialize_solution_step(self):
        super().initialize_solution_step()

        self.updated = False
        self.surrogate_updated = False

    def finalize_solution_step(self):
        super().finalize_solution_step()
        if not self.surrogate_updated:
            raise Exception('Surrogate not updated')
        if not self.updated:
            raise Exception('Not updated')

    def predict(self, x_in):
        x = x_in.copy()
        if self.updated:
            raise Exception('Already updated')
        if not self.surrogate_updated:
            raise Exception('Surrogate not updated: make sure that in the coupled solver'
                            'a surrogate is used that provides a surrogate solution')
        if self.predict_change:
            dxs = self.surrogate_data - self.surrogate_dataprev
            x.set_interface_data(self.dataprev + dxs)
        else:
            x.set_interface_data(self.surrogate_data)
        return x

    def update_surrogate(self, xs):
        if self.surrogate_updated:
            raise Exception('Surrogate already updated')
        self.surrogate_dataprev = self.surrogate_data
        self.surrogate_data = xs.get_interface_data()
        self.surrogate_updated = True

    def update(self, x):
        if self.updated:
            raise Exception('Already updated')
        self.dataprev = x.get_interface_data()
        self.updated = True

    def restart(self, restart_data):
        self.dataprev = restart_data['dataprev']
        self.surrogate_data = restart_data['surrogate_data']

    def check_restart_data(self, restart_data):
        model_type = restart_data['type']
        for key in ['predict_change']:
            new_value = self.settings.get(key)
            old_value = restart_data['settings'].get(key)
            if new_value != old_value:
                tools.print_info(f'"{model_type}" parameter "{key}" changed from {old_value} to {new_value}',
                                 layout='blue')

    def save_restart_data(self):
        return {'dataprev': self.dataprev, 'surrogate_data': self.surrogate_data}
