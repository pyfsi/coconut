from coconut import data_structure
from coconut.coupling_components.component import Component
from coconut.data_structure.interface import Interface
from coconut import tools

import os
from os.path import join
import numpy as np
import json



def create(parameters):
    return VolConsSolverV1(parameters)


class VolConsSolverV1(Component):
    @tools.time_initialize
    def __init__(self, parameters):
        super().__init__()

        self.env = None
        self.timestep = 0
        self.interface_list = None
        self.settings = parameters['settings']
        self.working_directory = join(os.getcwd(), self.settings['working_directory'])
        self.delta_t = self.settings['delta_t']
        self.timestep_start = self.settings['timestep_start']
        self.save_restart = self.settings.get('save_restart', 0)
        self.free_surface_updaters = []
        input_file_name = join(self.working_directory, self.settings['input_file'])
        with open(input_file_name, 'r') as parameter_file:
            input_parameters = json.load(parameter_file)

        #supply time parameters from coconut to the solver
        input_parameters['delta_t'] = self.delta_t
        input_parameters['time_step_start'] = self.timestep_start
        input_parameters['save_restart'] = self.save_restart
        # update the mesh path (rel to working directory)
        input_parameters['mesh_filename'] = os.path.join(self.working_directory, input_parameters['mesh_filename'])
        self.model = data_structure.Model()

        self.settings['interface_input']= []
        self.settings['interface_output']= []
        from vol_cons_solver.free_surface_updater import FreeSurfaceUpdater
        for mp_name in self.settings['interface_list']:
            input_parameters['sub_mesh_name'] = mp_name
            solver = FreeSurfaceUpdater(input_parameters)
            self.free_surface_updaters.append(solver)
            x0 = solver.get_mesh().points[:, 0]
            y0 = solver.get_mesh().points[:, 1]
            z0 = solver.get_mesh().points[:, 2]
            self.model.create_model_part(f'{mp_name}_input', x0, y0, z0, np.arange(0, x0.size))
            self.model.create_model_part(f'{mp_name}_output', x0, y0, z0, np.arange(0, x0.size))
            self.settings['interface_input'].append({'model_part': f'{mp_name}_input', 'variables': ['displacement']})
            self.settings['interface_output'].append({'model_part': f'{mp_name}_output', 'variables': ['pressure', 'traction']})

        # # Interfaces
        self.interface_input = Interface(self.settings['interface_input'], self.model)
        self.interface_output = Interface(self.settings['interface_output'], self.model)

        # time
        self.init_time = self.init_time
        self.run_time = 0.0

 
    def initialize(self):
        super().initialize()
        self.timestep = 0

    def initialize_solution_step(self):
        super().initialize_solution_step()
        for solver in self.free_surface_updaters:
            solver.initialize_time_step()
        self.timestep += 1

    @tools.time_solve_solution_step
    def solve_solution_step(self, interface_input):

        self.interface_input.set_interface_data(interface_input.get_interface_data())
        self.apply_input_displacement()
        self.update_output_pressure()
        return self.get_interface_output()

    def finalize_solution_step(self):
        super().finalize_solution_step()
        for solver in self.free_surface_updaters:
            solver.finalize_time_step()


    def finalize(self):
        super().finalize()


    def get_interface_input(self):
        return self.interface_input

    def get_interface_output(self):
        return self.interface_output

    def apply_input_displacement(self):
        for i, mp_name in enumerate(self.settings['interface_list']):
            input_mp_name = f'{mp_name}_input'
            displacement = self.interface_input.get_variable_data(input_mp_name, 'displacement')
            self.free_surface_updaters[i].apply_displacement(displacement)
            self.free_surface_updaters[i].update_free_surface()


    def update_output_pressure(self):

        for i, mp_name in enumerate(self.settings['interface_list']):
            output_mp_name = f'{mp_name}_output'
            pressure = self.free_surface_updaters[i].calculate_pressure()
            traction = np.zeros((pressure.size, 3))
            self.interface_output.set_variable_data(output_mp_name, 'pressure', pressure[:, np.newaxis])
            self.interface_output.set_variable_data(output_mp_name, 'traction', traction)

