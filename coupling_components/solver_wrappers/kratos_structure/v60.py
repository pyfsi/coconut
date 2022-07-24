from coconut.coupling_components.solver_wrappers.kratos_structure.base_solver_wrapper import \
    BaseSolverWrapperKratosStructure
from coconut import tools

import json
import os
import numpy as np
import re


def create(parameters):
    return SolverWrapperKratosStructure60(parameters)


class SolverWrapperKratosStructure60(BaseSolverWrapperKratosStructure):
    version = '60'

    def __init__(self, parameters):
        super().__init__(parameters)
        self.env = tools.get_solver_env(__name__, self.working_directory)

    def update_kratos_parameter_file(self, input_file_name):

        with open(input_file_name, 'r') as parameter_file:
            kratos_parameters = json.load(parameter_file)

        kratos_parameters['problem_data']['start_time'] = 0.0
        kratos_parameters['problem_data']['time_step'] = self.delta_t
        kratos_parameters['problem_data']['end_time'] = 1e15
        if 'structure_iterations' in self.settings:
            kratos_parameters['solver_settings']['max_iteration'] = self.settings['structure_iterations']
        kratos_parameters['interface_sub_model_parts_list'] = self.interface_sub_model_parts_list
        kratos_parameters['pressure_directions'] = self.check_pressure_directions()

        if self.save_restart:
            restart_save_dict = {'restart_control_type': 'step',
                                 'restart_save_frequency': abs(self.save_restart),
                                 # kratos 6.0 does not support negative numbers
                                 'save_restart': True}
            kratos_parameters['solver_settings']['restart_settings'] = restart_save_dict

        if self.timestep_start != 0:
            restart_load_dict = {'load_restart': True,
                                 'restart_load_file_label': f'{self.timestep_start}'}
            kratos_parameters['solver_settings']['restart_settings'].update(restart_load_dict)

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
                        search_string = r'\n' + variable + r' CRITERION.*[nN]orm = ' + r'(' + float_pattern + r')'
                        var_residual_list = re.findall(search_string, iteration_block)
                        if var_residual_list:
                            # last initial residual of the non-linear iteration
                            var_residual = float(var_residual_list[-1])
                            residual_array[i] = var_residual
                        else:
                            raise RuntimeError(f'{variable} CRITERION not found in kratos log file')

                    with open(self.res_filepath, 'a') as f:
                        np.savetxt(f, [residual_array], delimiter=', ')
