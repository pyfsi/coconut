from coconut.coupling_components.solver_wrappers.kratos_structure.kratos_structure import SolverWrapperKratosStructure
from coconut import tools

import json
import os
import numpy as np
import re


def create(parameters):
    return SolverWrapperKratosStructure91(parameters)


class SolverWrapperKratosStructure91(SolverWrapperKratosStructure):
    version = '91'

    def __init__(self, parameters):
        super().__init__(parameters)
        self.env = tools.get_solver_env(__name__, self.working_directory)
        self.check_software()

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
                            pass
			    #raise RuntimeError(f'{variable} or {variable} CRITERION not found in kratos log file')

                    with open(self.res_filepath, 'a') as f:
                        np.savetxt(f, [residual_array], delimiter=', ')
