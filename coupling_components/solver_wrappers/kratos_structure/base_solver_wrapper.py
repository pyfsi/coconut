import os
import shutil
import time
from os.path import join
from subprocess import Popen

import numpy as np
import pandas as pd
from coconut import data_structure
from coconut import tools
from coconut.coupling_components.component import Component
from coconut.data_structure.interface import Interface


def create(parameters):
    return BaseSolverWrapperKratosStructure(parameters)


class BaseSolverWrapperKratosStructure(Component):
    # version specific parameters
    version = None  # KratosMultiphysics version, set in sub-class, for version 9.1 f. ex.: '91'

    @tools.time_initialize
    def __init__(self, parameters):
        super().__init__()

        if self.version is None:
            raise NotImplementedError(
                'Base class method called, class variable version needs to be set in the derived class')

        self.settings = parameters['settings']
        self.working_directory = join(os.getcwd(), self.settings['working_directory'])
        self.env = None

        self.delta_t = self.settings['delta_t']
        self.timestep_start = self.settings['timestep_start']
        self.timestep = self.timestep_start
        self.save_restart = self.settings.get('save_restart', 0)
        self.iteration = None

        self.interface_sub_model_parts_list = self.settings['kratos_interface_sub_model_parts_list']

        self.model = None
        self.interface_input = None
        self.interface_output = None
        self.kratos_process = None

        self.check_interface()

        # residuals
        self.residual_variables = self.settings.get('residual_variables', None)
        self.res_filepath = os.path.join(self.working_directory, 'residuals.csv')

        if self.residual_variables is not None:
            self.write_residuals_fileheader()

        # coupling convergence
        self.coupling_convergence = False

        # time
        self.init_time = self.init_time
        self.run_time = 0.0

        # debug
        self.debug = self.settings.get('debug', False)  # save copy of input and output files in every iteration

    @tools.time_initialize
    def initialize(self):
        super().initialize()

        input_file_name = join(self.working_directory, self.settings['input_file'])
        self.update_kratos_parameter_file(input_file_name)

        self.model = data_structure.Model()

        dir_path = os.path.dirname(os.path.realpath(__file__))
        run_script_file = os.path.join(dir_path, f'run_kratos_structural_{self.version}.py')

        self.kratos_process = Popen(f'python3 {run_script_file} {input_file_name} &> log',
                                    shell=True, cwd=self.working_directory, env=self.env)
        self.wait_message('start_ready')

        for mp_name in self.interface_sub_model_parts_list:
            file_path = os.path.join(self.working_directory, f'{mp_name}_nodes.csv')
            node_data = pd.read_csv(file_path, skipinitialspace=True)
            node_ids = np.array(node_data.node_id)
            x0 = np.array(node_data.x0)
            y0 = np.array(node_data.y0)
            z0 = np.array(node_data.z0)
            self.model.create_model_part(f'{mp_name}_input', x0, y0, z0, node_ids)
            self.model.create_model_part(f'{mp_name}_output', x0, y0, z0, node_ids)

        # interfaces
        self.interface_input = Interface(self.settings['interface_input'], self.model)
        self.interface_output = Interface(self.settings['interface_output'], self.model)

    def initialize_solution_step(self):
        super().initialize_solution_step()

        self.iteration = 0
        self.timestep += 1

        self.send_message('next')
        self.wait_message('next_ready')

    @tools.time_solve_solution_step
    def solve_solution_step(self, interface_input):
        self.iteration += 1

        self.interface_input.set_interface_data(interface_input.get_interface_data())
        self.write_input_data()
        self.send_message('continue')
        self.wait_message('continue_ready')
        self.update_interface_output()
        return self.get_interface_output()

    def finalize_solution_step(self):
        super().finalize_solution_step()
        self.send_message('save')
        self.wait_message('save_ready')
        if self.residual_variables is not None:
            self.write_residuals()

    def finalize(self):
        super().finalize()
        self.send_message('stop')
        self.wait_message('stop_ready')
        self.remove_all_messages()
        self.kratos_process.wait()

    def get_interface_input(self):
        return self.interface_input

    def get_interface_output(self):
        return self.interface_output

    def write_input_data(self):
        interface_sub_model_parts_list = self.settings['kratos_interface_sub_model_parts_list']

        for mp_name in interface_sub_model_parts_list:
            input_mp_name = f'{mp_name}_input'
            input_mp = self.model.get_model_part(input_mp_name)
            node_ids = np.array([input_mp.id[i] for i in range(input_mp.size)])

            file_path_pr = os.path.join(self.working_directory, f'{mp_name}_pressure.csv')
            pressure_array = self.interface_input.get_variable_data(input_mp_name, 'pressure')
            pressure_df = pd.DataFrame({'node_id': node_ids, 'pressure': pressure_array[:, 0]})
            pressure_df.to_csv(file_path_pr, index=False)

            file_path_tr = os.path.join(self.working_directory, f'{mp_name}_traction.csv')
            traction_array = self.interface_input.get_variable_data(input_mp_name, 'traction')
            traction_df = pd.DataFrame({'node_id': node_ids, 'traction_x': traction_array[:, 0],
                                        'traction_y': traction_array[:, 1], 'traction_z': traction_array[:, 2]})
            traction_df.to_csv(file_path_tr, index=False)

            # copy input data for debugging
            if self.debug:
                pressure_df.to_csv(file_path_pr[:-4] + f'_ts{self.timestep}_it{self.iteration}.csv', index=False)
                traction_df.to_csv(file_path_tr[:-4] + f'_ts{self.timestep}_it{self.iteration}.csv', index=False)

    def update_interface_output(self):
        interface_sub_model_parts_list = self.settings['kratos_interface_sub_model_parts_list']

        for mp_name in interface_sub_model_parts_list:
            output_mp_name = f'{mp_name}_output'
            file_path = os.path.join(self.working_directory, f'{mp_name}_displacement.csv')
            disp_data = pd.read_csv(file_path, skipinitialspace=True)
            disp_x = np.array(disp_data.displacement_x)
            disp_y = np.array(disp_data.displacement_y)
            disp_z = np.array(disp_data.displacement_z)
            displacement = np.column_stack((disp_x, disp_y, disp_z))
            self.interface_output.set_variable_data(output_mp_name, 'displacement', displacement)

            if self.debug:
                shutil.copy(file_path, file_path[:-4] + f'_ts{self.timestep}_it{self.iteration}.csv')

    def update_kratos_parameter_file(self, input_file_name):
        raise NotImplementedError('Base class method is called, should be implemented in sub-class')

    def check_interface(self):
        input_interface_model_parts = [param['model_part'] for param in self.settings['interface_input']]
        output_interface_model_parts = [param['model_part'] for param in self.settings['interface_output']]
        sub_mp_name_list = self.settings['kratos_interface_sub_model_parts_list']

        for sub_mp_name in sub_mp_name_list:
            if f'{sub_mp_name}_input' not in input_interface_model_parts:
                raise RuntimeError(
                    f'Error in json file: {sub_mp_name}_input not listed in "interface_input": '
                    f'{self.settings["interface_input"]}.\n. <sub_mp_name> in the '
                    f'"kratos_interface_sub_model_parts_list" in json file should have corresponding '
                    f'<sub_mp_name>_input in "interface_input" list. ')

            if f'{sub_mp_name}_output' not in output_interface_model_parts:
                raise RuntimeError(
                    f'Error in json file: {sub_mp_name}_output not listed in "interface_output": '
                    f'{self.settings["interface_output"]}.\n. <sub_mp_name> in the '
                    f'"kratos_interface_sub_model_parts_list" in json file should have corresponding '
                    f'<sub_mp_name>_output in "interface_output" list.')

    def check_pressure_directions(self):
        pressure_directions = self.settings.get('pressure_directions')

        if pressure_directions is None:
            return [1] * len(self.interface_sub_model_parts_list)  # default pressure in direction of normal
        if len(pressure_directions) != len(self.interface_sub_model_parts_list):
            raise ValueError(f'Error in json file: "pressure_directions" {pressure_directions} should have length '
                             f'{len(self.interface_sub_model_parts_list)}')
        if [abs(i) for i in pressure_directions] != [1] * len(self.interface_sub_model_parts_list):
            raise ValueError(
                f'Error in json file: "pressure_directions" {pressure_directions} must only contain 1 or -1')
        return pressure_directions

    def send_message(self, message):
        file = join(self.working_directory, message + '.coco')
        open(file, 'w').close()
        return

    def wait_message(self, message):
        file = join(self.working_directory, message + '.coco')
        while not os.path.isfile(file):
            time.sleep(0.01)
        os.remove(file)
        return

    def check_message(self, message):
        file = join(self.working_directory, message + '.coco')
        if os.path.isfile(file):
            os.remove(file)
            return True
        return False

    def remove_all_messages(self):
        for file_name in os.listdir(self.working_directory):
            if file_name.endswith('.coco'):
                file = join(self.working_directory, file_name)
                os.remove(file)

    def write_residuals_fileheader(self):
        header = ''
        sep = ', '
        with open(self.res_filepath, 'w') as f:
            f.write('# Residuals\n')
            for variable in self.residual_variables:
                header += variable + sep
            f.write(header.strip(sep) + '\n')

    def write_residuals(self):
        raise NotImplementedError('Base class method is called, should be implemented in sub-class')
