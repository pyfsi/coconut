from coconut.tools import create_instance
from coconut.coupling_components.component import Component
from coconut import tools
from coconut.data_structure import Interface


def create(parameters):
    return SolverWrapperBaffle(parameters)


class SolverWrapperBaffle(Component):
    def __init__(self, parameters):
        super().__init__()

        self.parameters = parameters
        self.settings = parameters['settings']
        self.save_restart = self.settings.get('save_restart', 0)  # time step interval to save restart data
        self.settings['save_restart'] = self.save_restart

        # create solvers
        sol_wrapper_param = self.settings['solver_wrapper']
        tools.pass_on_parameters(self.settings, sol_wrapper_param['settings'], ['timestep_start', 'delta_t', 'save_restart'])
        self.solver_wrapper = create_instance(sol_wrapper_param)

        mapper_settings = {"type": "mappers.nearest",
                           "settings": {
                               "directions": ["x", "y", "z"]
                           }}

        # create mappers
        self.mapper_input = create_instance(mapper_settings)
        self.mapper_output = create_instance(mapper_settings)
        self.master_mp_inp_name = self.settings['master_model_part_inp']
        self.master_mp_out_name = self.settings['master_model_part_out']
        self.slave_mp_inp_name = self.settings['slave_model_part_inp']
        self.slave_mp_out_name = self.settings['slave_model_part_out']
        self.interface_input = None
        self.interface_output = None

        # time
        self.init_time = 0.0
        self.run_time = 0.0

    def initialize(self):
        super().initialize()
        self.solver_wrapper.initialize()
        interface_input = self.solver_wrapper.get_interface_input().copy()
        interface_output = self.solver_wrapper.get_interface_output().copy()
        interface_input_param = interface_input.parameters
        interface_output_param = interface_output.parameters
        for param in interface_input_param:
            if self.slave_mp_inp_name in param['model_part']:
                interface_input_param.remove(param)
        for param in interface_output_param:
            if self.slave_mp_out_name in param['model_part']:
                interface_output_param.remove(param)

        self.interface_input = Interface(interface_input_param, interface_input.model)
        self.interface_output = Interface(interface_output_param, interface_output.model)
        mp_inp_from = self.interface_input.get_model_part(self.master_mp_inp_name)
        mp_inp_to = interface_input.get_model_part(self.slave_mp_inp_name)
        mp_out_from = interface_output.get_model_part(self.slave_mp_out_name)
        mp_out_to = self.interface_output.get_model_part(self.master_mp_out_name)
        self.mapper_input.initialize(mp_inp_from, mp_inp_to)
        self.mapper_output.initialize(mp_out_from, mp_out_to)

    def initialize_solution_step(self):
        super().initialize_solution_step()
        self.solver_wrapper.initialize_solution_step()

    @tools.time_solve_solution_step
    def solve_solution_step(self, interface_input):
        sol_wra_inp_int = self.apply_input_data_to_sol_wrapper(interface_input, 'displacement')
        sol_wra_out_int = self.solver_wrapper.solve_solution_step(sol_wra_inp_int)
        self.apply_output_data_from_sol_wrapper(sol_wra_out_int, 'pressure', -1)
        self.apply_output_data_from_sol_wrapper(sol_wra_out_int, 'traction', 1)
        return self.interface_output

    def finalize_solution_step(self):
        super().finalize_solution_step()
        self.solver_wrapper.finalize_solution_step()

    def finalize(self):
        super().finalize()
        self.solver_wrapper.finalize()

    def output_solution_step(self):
        super().output_solution_step()
        self.solver_wrapper.output_solution_step()

    def get_interface_input(self):
        # does not contain most recent data
        return self.interface_input

    def get_interface_output(self):
        # does not contain most recent data
        return self.interface_output

    def apply_input_data_to_sol_wrapper(self, interface_input, variable):
        sol_wra_inp_int =  self.solver_wrapper.get_interface_input()
        master_inp_var_data = interface_input.get_variable_data(self.master_mp_inp_name, variable)
        sol_wra_inp_int.set_variable_data(self.master_mp_inp_name, variable, master_inp_var_data)
        master_inp_args_from = (interface_input, self.master_mp_inp_name, variable)
        slave_inp_args_to = (sol_wra_inp_int, self.slave_mp_inp_name, variable)
        self.mapper_input(master_inp_args_from, slave_inp_args_to)
        return sol_wra_inp_int

    def apply_output_data_from_sol_wrapper(self, sol_wrapper_interface_output, variable, mul_factor):
        master_out_var_data = sol_wrapper_interface_output.get_variable_data(self.master_mp_out_name, variable)
        self.interface_output.set_variable_data(self.master_mp_out_name, variable, master_out_var_data)
        interface_output_to = self.interface_output.copy()
        interface_output_to *= 0
        slave_out_args = (sol_wrapper_interface_output, self.slave_mp_out_name, variable)
        master_out_args = (interface_output_to,self.master_mp_out_name, variable)
        self.mapper_output(slave_out_args, master_out_args)
        self.interface_output += mul_factor*interface_output_to

    def print_components_info(self, pre):
        tools.print_info(pre, 'The component ', self.__class__.__name__, ' combines the interface with baffles for a solver wrapper:')
        pre = tools.update_pre(pre)
        tools.print_info(pre, '└─', 'solver wrapper:')
        self.solver_wrapper.print_components_info(pre + '  └─')
