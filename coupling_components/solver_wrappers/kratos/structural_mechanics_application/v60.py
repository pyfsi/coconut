from coconut import data_structure
from coconut.coupling_components.component import Component
from coconut.data_structure.interface import Interface
from coconut import tools

import json
import os
import time
from os.path import join
from subprocess import Popen
import pandas as pd
import numpy as np
import re


def create(parameters):
    return SolverWrapperKratosStructure60(parameters)


class SolverWrapperKratosStructure60(Component):
    def __init__(self, parameters):
        super().__init__()

        self.settings = parameters["settings"]
        self.working_directory = join(os.getcwd(), self.settings["working_directory"])
        self.env = tools.get_solver_env(__name__, self.working_directory)
        delta_t = self.settings["delta_t"]
        timestep_start = self.settings["timestep_start"]
        dimensions = self.settings["dimensions"]

        input_file_name = join(self.working_directory, self.settings["input_file"])

        with open(input_file_name, "r") as parameter_file:
            kratos_parameters = json.load(parameter_file)

        kratos_parameters["problem_data"]["start_time"] = timestep_start
        kratos_parameters["problem_data"]["time_step"] = delta_t
        kratos_parameters["problem_data"]["domain_size"] = dimensions
        kratos_parameters["problem_data"]["end_time"] = 1e15

        interface_sub_model_parts_list = self.settings["kratos_interface_sub_model_parts_list"]

        kratos_parameters["interface_sub_model_parts_list"] = interface_sub_model_parts_list

        with open(os.path.join(self.working_directory, input_file_name), 'w') as f:
            json.dump(kratos_parameters, f, indent=4)

        self.check_interface()

        self.model = data_structure.Model()

        dir_path = os.path.dirname(os.path.realpath(__file__))
        run_script_file = os.path.join(dir_path, 'run_kratos_structural_60.py')

        self.kratos_process = Popen(f'python3 {run_script_file} {input_file_name} &> log',
                                    shell=True, cwd=self.working_directory, env=self.env)

        self.wait_message('start_ready')

        for mp_name in interface_sub_model_parts_list:
            file_path = os.path.join(self.working_directory, f'{mp_name}_nodes.csv')
            node_data = pd.read_csv(file_path, skipinitialspace=True)
            node_ids = np.array(node_data.node_id)
            x0 = np.array(node_data.x0)
            y0 = np.array(node_data.y0)
            z0 = np.array(node_data.z0)
            self.model.create_model_part(f'{mp_name}_input', x0, y0, z0, node_ids)
            self.model.create_model_part(f'{mp_name}_output', x0, y0, z0, node_ids)

        # # Interfaces
        self.interface_input = Interface(self.settings["interface_input"], self.model)
        self.interface_output = Interface(self.settings["interface_output"], self.model)

        # run time
        self.run_time = 0.0

        self.residual_variables = self.settings.get('residual_variables', None)
        self.res_filepath = os.path.join(self.working_directory, 'residuals.dat')

        if not self.residual_variables is None:
            self.write_residuals_fileheader()

    def initialize(self):
        super().initialize()
        self.timestep = 0

    def initialize_solution_step(self):
        super().initialize_solution_step()
        self.timestep +=1

        self.send_message('next')
        self.wait_message('next_ready')

    @tools.time_solve_solution_step
    def solve_solution_step(self, interface_input):

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
        if not self.residual_variables is None:
            self.write_residuals()

    def finalize(self):
        super().finalize()
        self.send_message('stop')
        self.wait_message('stop_ready')
        self.remove_all_messages()
        self.kratos_process.kill()

    def get_interface_input(self):
        return self.interface_input

    def get_interface_output(self):
        return self.interface_output

    def write_input_data(self):
        interface_sub_model_parts_list = self.settings["kratos_interface_sub_model_parts_list"]

        for mp_name in interface_sub_model_parts_list:
            input_mp_name = f'{mp_name}_input'
            input_mp = self.model.get_model_part(input_mp_name)
            file_path_pr = os.path.join(self.working_directory, f'{mp_name}_pressure.csv')
            with open(file_path_pr, 'w') as f:
                f.write('node_id, pressure\n')
            file_path_sl = os.path.join(self.working_directory, f'{mp_name}_surface_load.csv')
            with open(file_path_sl, 'w') as f:
                f.write('node_id, surface_load_x, surface_load_y, surface_load_z\n')

            pressure_array = np.ravel(self.interface_input.get_variable_data(input_mp_name, 'pressure'))
            surface_load_array = self.interface_input.get_variable_data(input_mp_name, 'traction')

            for i in range(0, input_mp.size):
                with open(file_path_pr, 'a') as f:
                    f.write(str(input_mp.id[i]) + ', ' + str(pressure_array[i]) + '\n')

                with open(file_path_sl, 'a') as f:
                    f.write(str(input_mp.id[i]) + ', ' + str(surface_load_array[i, 0]) + ', ' + str(
                        surface_load_array[i, 1]) + ', ' + str(
                        surface_load_array[i, 2]) + '\n')

    def update_interface_output(self):
        interface_sub_model_parts_list = self.settings["kratos_interface_sub_model_parts_list"]

        for mp_name in interface_sub_model_parts_list:
            output_mp_name = f'{mp_name}_output'
            file_path = os.path.join(self.working_directory, f'{mp_name}_displacement.csv')
            disp_data = pd.read_csv(file_path, skipinitialspace=True)
            disp_x = np.array(disp_data.displacement_x)
            disp_y = np.array(disp_data.displacement_y)
            disp_z = np.array(disp_data.displacement_z)
            displacement = np.column_stack((disp_x, disp_y, disp_z))
            self.interface_output.set_variable_data(output_mp_name, 'displacement', displacement)

    def check_interface(self):

        input_interface_model_parts = [param["model_part"] for param in self.settings["interface_input"]]
        output_interface_model_parts = [param["model_part"] for param in self.settings["interface_output"]]
        sub_mp_name_list = self.settings["kratos_interface_sub_model_parts_list"]

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

    def send_message(self, message):
        file = join(self.working_directory, message + ".coco")
        open(file, 'w').close()
        return

    def wait_message(self, message):
        file = join(self.working_directory, message + ".coco")
        while not os.path.isfile(file):
            time.sleep(0.01)
        os.remove(file)
        return

    def check_message(self, message):
        file = join(self.working_directory, message + ".coco")
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
        header = '# '
        nr_spaces =15
        for var in self.residual_variables:
            header += var + ' '*nr_spaces
        with open(self.res_filepath, 'w') as f:
            f.write('# Residuals\n')
            f.write(header + '\n')

    def write_residuals(self):
        nr_spaces = 15
        float_pattern = r'[+-]?\d*\.?\d*[eE]?[+-]?\d*'
        log_filepath = os.path.join(self.working_directory, f'log')
        if os.path.isfile(log_filepath):
            with open(log_filepath, 'r') as f:
                log_string = f.read()
            time_start_string = r'STEP:\s+' + str(self.timestep -1)
            time_end_string = r'STEP:\s+' + str(self.timestep)
            match = re.search(time_start_string + r'(.*)' + time_end_string, log_string, flags=re.S)
            if not match is None:
                time_block = match.group(1)
                iteration_block_list = re.findall(
                    r'Coupling iteration: \d+' + r'(.*?)' + r'Coupling iteration \d+ end', time_block, flags=re.S)
                for iteration_block in iteration_block_list:
                    for variable in self.residual_variables:
                        search_string = r'\n' + variable + r' CRITERION.*[Nn]orm = ' + r'(' + float_pattern + r')'
                        var_residual_list = re.findall(search_string, iteration_block)
                        if var_residual_list:
                            # last initial residual of the non-linear iteration
                            var_residual = float(var_residual_list[-1])
                        else:
                            raise RuntimeError(f'{variable} CRITERION not found in kratos log file')

                        with open(self.res_filepath, 'a') as f:
                            f.write(str(var_residual) + ' ' * nr_spaces)
                    with open(self.res_filepath, 'a') as f:
                        f.write('\n')
