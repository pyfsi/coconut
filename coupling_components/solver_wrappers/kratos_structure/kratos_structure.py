from coconut.coupling_components.solver_wrappers.solver_wrapper import SolverWrapper
from coconut.data_structure.interface import Interface, Model
from coconut import tools

import shutil
import os
from os.path import join
from subprocess import Popen, run, PIPE
import pandas as pd
import numpy as np
import json
import re


def create(parameters):
    return SolverWrapperKratosStructure(parameters)


class SolverWrapperKratosStructure(SolverWrapper):
    version = None  # KratosMultiphysics version, set in subclass, for version 9.1 f. ex.: '91'
    check_coupling_convergence_possible = True  # can solver check convergence after 1 iteration?

    # define input and output variables
    accepted_in_var = ['pressure', 'traction']
    accepted_out_var = ['displacement']

    @tools.time_initialize
    def __init__(self, parameters):
        super().__init__(parameters)

        if self.version is None:
            raise NotImplementedError(
                'Base class method called, class variable version needs to be set in the derived class')

        self.settings = parameters['settings']
        self.working_directory = join(os.getcwd(), self.settings['working_directory'])
        self.env = None
        self.coco_messages = tools.CocoMessages(self.working_directory)
        self.coco_messages.remove_all_messages()

        self.delta_t = self.settings['delta_t']
        self.timestep_start = self.settings['timestep_start']
        self.timestep = self.timestep_start
        self.save_restart = self.settings.get('save_restart', 0)
        self.iteration = None

        self.interface_sub_model_parts_list = self.settings['kratos_interface_sub_model_parts_list']

        self.model = None
        self.kratos_process = None

        self.check_interface()

        # residuals
        self.residual_variables = self.settings.get('residual_variables', None)
        self.res_filepath = os.path.join(self.working_directory, 'residuals.csv')

        if self.residual_variables is not None:
            self.write_residuals_fileheader()

    @tools.time_initialize
    def initialize(self):
        super().initialize()

        input_file_name = join(self.working_directory, self.settings['input_file'])
        self.update_kratos_parameter_file(input_file_name)

        self.model = Model()

        dir_path = os.path.dirname(os.path.realpath(__file__))
        run_script_file = os.path.join(dir_path, f'run_kratos_structural_v{self.version}.py')

        self.kratos_process = Popen(f'python3 {run_script_file} {input_file_name} &> kratos.log',
                                    shell=True, cwd=self.working_directory, env=self.env)

        # pass on process to coco_messages for polling
        self.coco_messages.set_process(self.kratos_process)

        self.coco_messages.wait_message('start_ready')

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

        self.coco_messages.send_message('next')
        self.coco_messages.wait_message('next_ready')

    @tools.time_solve_solution_step
    def solve_solution_step(self, interface_input):
        self.iteration += 1

        self.interface_input.set_interface_data(interface_input.get_interface_data())
        self.write_input_data()
        self.coco_messages.send_message('continue')
        self.coco_messages.wait_message('continue_ready')

        if self.check_coupling_convergence:
            # check if Kratos converged after one iteration
            self.coupling_convergence = self.coco_messages.check_message('solver_converged')
            if self.print_coupling_convergence and self.coupling_convergence:
                tools.print_info(f'{self.__class__.__name__} converged')

        self.update_interface_output()
        return self.get_interface_output()

    def finalize_solution_step(self):
        super().finalize_solution_step()

    @tools.time_save
    def output_solution_step(self):
        super().output_solution_step()

        self.coco_messages.send_message('save')
        self.coco_messages.wait_message('save_ready')
        if self.residual_variables is not None:
            self.write_residuals()

    def finalize(self):
        super().finalize()

        self.coco_messages.send_message('stop')
        self.coco_messages.wait_message('stop_ready')
        self.coco_messages.remove_all_messages()
        self.kratos_process.wait()

        # remove unnecessary files
        with tools.cd(self.working_directory):
            for mp_name in self.interface_sub_model_parts_list:
                os.remove(f'{mp_name}_displacement.csv')
                os.remove(f'{mp_name}_pressure.csv')
                os.remove(f'{mp_name}_traction.csv')
                os.remove(f'{mp_name}_nodes.csv')

    def write_input_data(self):
        for mp_name in self.interface_sub_model_parts_list:
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
        for mp_name in self.interface_sub_model_parts_list:
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

        with open(input_file_name, 'r') as parameter_file:
            kratos_parameters = json.load(parameter_file)

        kratos_parameters['problem_data']['start_time'] = 0.0
        kratos_parameters['solver_settings']['time_stepping']['time_step'] = self.delta_t
        kratos_parameters['problem_data']['end_time'] = 1e15
        if 'structure_iterations' in self.settings:
            kratos_parameters['solver_settings']['max_iteration'] = self.settings['structure_iterations']
        kratos_parameters['interface_sub_model_parts_list'] = self.interface_sub_model_parts_list
        kratos_parameters['pressure_directions'] = self.check_pressure_directions()
        kratos_parameters['check_coupling_convergence'] = self.check_coupling_convergence

        if self.save_restart:
            restart_save_dict = {'restart_processes': [{'python_module': 'save_restart_process',
                                                        'kratos_module': 'KratosMultiphysics',
                                                        'process_name': 'SaveRestartProcess',
                                                        'Parameters': {
                                                            'model_part_name': 'Structure',
                                                            'restart_control_type': 'step',
                                                            'restart_save_frequency': abs(self.save_restart)}}]}
            if self.save_restart < 0:
                restart_save_dict['restart_processes'][0]['Parameters']['max_files_to_keep'] = 1
            kratos_parameters['output_processes'].update(restart_save_dict)

        if self.timestep_start != 0:
            restart_load_dict = {'restart_load_file_label': str(self.timestep_start),
                                 'input_type': 'rest',
                                 'input_filename': 'Structure'}
            kratos_parameters['solver_settings']['model_import_settings'].update(restart_load_dict)

        with open(os.path.join(self.working_directory, input_file_name), 'w') as f:
            json.dump(kratos_parameters, f, indent=2)

    def check_software(self):
        with open('check_software.py', 'w') as f:
            f.write('import KratosMultiphysics\nfrom KratosMultiphysics import StructuralMechanicsApplication')
        process = run('python3 check_software.py', stdout=PIPE, stderr=PIPE, encoding='utf-8', shell=True, env=self.env)
        os.remove('check_software.py')
        if process.returncode != 0:
            raise RuntimeError(f'KratosMultiphysics not loaded properly. Check if the solver load commands for '
                               f'the "machine_name" are correct in solver_modules.py.')

        # check version
        line = process.stdout.splitlines()[4]   # fifth line contains 'Multi-Physics X.X' with XX the version number
        stripped_version = line.split('Multi-Physics ')[1]
        version = stripped_version[0:3:2]
        if version != self.version:
            raise RuntimeError(f'Kratos {self.version[0]}.{self.version[1]} should be loaded!'
                               f' Currently, Kratos {stripped_version} is loaded. Check if the solver load commands'
                               f' for the "machine_name" are correct in solver_modules.py.')

    def check_interface(self):
        input_interface_model_parts = [param['model_part'] for param in self.settings['interface_input']]
        output_interface_model_parts = [param['model_part'] for param in self.settings['interface_output']]

        for sub_mp_name in self.interface_sub_model_parts_list:
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
        pressure_directions = self.settings['pressure_directions']

        if len(pressure_directions) != len(self.interface_sub_model_parts_list):
            raise ValueError(f'Error in json file: "pressure_directions" {pressure_directions} should have length '
                             f'{len(self.interface_sub_model_parts_list)}')
        if [abs(i) for i in pressure_directions] != [1] * len(self.interface_sub_model_parts_list):
            raise ValueError(
                f'Error in json file: "pressure_directions" {pressure_directions} must only contain 1 or -1')
        return pressure_directions

    def write_residuals_fileheader(self):
        header = ''
        sep = ', '
        with open(self.res_filepath, 'w') as f:
            f.write('# Residuals\n')
            for variable in self.residual_variables:
                header += variable + sep
            f.write(header.strip(sep) + '\n')

    def write_residuals(self):
        float_pattern = r'[+-]?\d*\.?\d*[eE]?[+-]?\d*'
        log_filepath = os.path.join(self.working_directory, f'log')
        if os.path.isfile(log_filepath):
            with open(log_filepath, 'r') as f:
                log_string = f.read()
            time_start_string = r'STEP:\s+' + str(self.timestep - 1)
            time_end_string = r'STEP:\s+' + str(self.timestep)
            match = re.search(time_start_string + r'(.*)' + time_end_string, log_string, flags=re.S)
            if match is not None:
                time_block = match.group(1)
                iteration_block_list = re.findall(
                    r'Coupling iteration: \d(.*?)Coupling iteration \d+ end', time_block, flags=re.S)
                for iteration_block in iteration_block_list:
                    residual_array = np.empty(len(self.residual_variables))
                    for i, variable in enumerate(self.residual_variables):
                        search_string = r'\n.*' + variable + r' CRITERION.*[nN]orm = +' + r'(' + float_pattern + r')'
                        var_residual_list_1 = re.findall(search_string, iteration_block)
                        search_string = r'\n ' + variable + r'.*abs = ' + r'(' + float_pattern + r')'
                        var_residual_list_2 = re.findall(search_string, iteration_block)
                        if var_residual_list_1:
                            # last initial residual of the non-linear iteration
                            var_residual = float(var_residual_list_1[-1])
                            residual_array[i] = var_residual
                        elif var_residual_list_2:
                            # last initial residual of the non-linear iteration
                            var_residual = float(var_residual_list_2[-1])
                            residual_array[i] = var_residual
                        else:
                            raise RuntimeError(f'{variable} or {variable} CRITERION not found in kratos log file')

                    with open(self.res_filepath, 'a') as f:
                        np.savetxt(f, [residual_array], delimiter=', ')
