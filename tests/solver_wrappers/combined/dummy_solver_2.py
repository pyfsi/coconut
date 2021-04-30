from coconut import data_structure
from coconut.coupling_components.component import Component
from coconut.data_structure.interface import Interface
from coconut import tools

import numpy as np

lx = 1.0
ly = 1.0
nx = 11
ny = 11


def create(parameters):
    return DummySolver2(parameters)


class DummySolver2(Component):
    def __init__(self, parameters):
        super().__init__()
        self.settings = parameters['settings']
        self.model = data_structure.Model()
        dx = lx / (nx - 1)
        dy = ly / (ny - 1)
        perturb_factor = 0.0
        x = np.linspace(0, lx, nx) + np.random.rand(nx) * dx * perturb_factor
        y = np.linspace(0, ly, ny) + np.random.rand(nx) * dy * perturb_factor

        xx, yy = np.meshgrid(x, y)
        x0_in = xx.ravel()
        y0_in = yy.ravel()

        x0_out = xx.ravel()
        y0_out = yy.ravel()

        z0_in = np.zeros_like(x0_in)
        z0_out = np.zeros_like(x0_out)

        node_ids_in = np.arange(0, nx * ny)
        node_ids_out = np.arange(0, nx * ny)

        self.mp_name_in_list = []
        for interface_settings in self.settings['interface_input']:
            self.model.create_model_part(interface_settings['model_part'], x0_in, y0_in, z0_in, node_ids_in)
            self.mp_name_in_list.append(interface_settings['model_part'])
        self.mp_name_out_list = []
        for interface_settings in self.settings['interface_output']:
            self.model.create_model_part(interface_settings['model_part'], x0_out, y0_out, z0_out, node_ids_out)
            self.mp_name_out_list.append(interface_settings['model_part'])

        # # Interfaces
        self.interface_input = Interface(self.settings["interface_input"], self.model)
        self.interface_output = Interface(self.settings["interface_output"], self.model)

        # run time
        self.run_time = 0.0

    def initialize(self):
        super().initialize()
        self.timestep = 0

    def initialize_solution_step(self):
        super().initialize_solution_step()
        self.timestep += 1

    @tools.time_solve_solution_step
    def solve_solution_step(self, interface_input):
        interface_data = interface_input.get_interface_data()
        for mp_name in self.mp_name_out_list:
            pressure_array, traction_array = self.calculate_output(interface_data, self.interface_output, mp_name)
            self.interface_output.set_variable_data(mp_name, 'pressure', pressure_array)
            self.interface_output.set_variable_data(mp_name, 'traction', traction_array)

        return self.get_interface_output()

    def finalize_solution_step(self):
        super().finalize_solution_step()

    def finalize(self):
        super().finalize()

    def get_interface_input(self):
        return self.interface_input

    def get_interface_output(self):
        return self.interface_output

    def calculate_output(self, data, out_interface, out_mp_name):
        nr_nodes = out_interface.get_model_part(out_mp_name).size
        norm = np.linalg.norm(data)
        min = np.min(data)
        pressure = norm
        traction = [min, -1 * min, 2 * min]
        return np.full((nr_nodes, 1), pressure), np.full((nr_nodes, 3), traction)


if __name__ == '__main__':
    import json

    with open('parameters.json', 'r') as parameter_file:
        parameters = json.load(parameter_file)

    sol_param = parameters['settings']['solver_wrappers'][1]
    dummy_solver = create(sol_param)
    print(dummy_solver.get_interface_output().get_interface_data())
