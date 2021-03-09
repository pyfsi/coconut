from coconut.tools import create_instance
from coconut.coupling_components.component import Component
from coconut import tools
from coconut import data_structure

import numpy as np


def create(parameters):
    return SolverWrapperMapped_update(parameters)


class SolverWrapperMapped_update(Component):
    def __init__(self, parameters):
        super().__init__()

        self.parameters = parameters
        self.settings = parameters["settings"]

        # create solver
        tools.pass_on_parameters(self.settings, self.settings["solver_wrapper"]["settings"],
                                 ["timestep_start", "delta_t"])
        self.solver_wrapper = create_instance(self.settings["solver_wrapper"])

        # create mappers
        self.mapper_interface_input = create_instance(self.settings["mapper_interface_input"])
        self.mapper_interface_output = create_instance(self.settings["mapper_interface_output"])

        self.mapper_interface_input_update_to = create_instance(self.settings["mapper_interface_input_update_to"])
        self.mapper_interface_input_update_from = create_instance(self.settings["mapper_interface_input_update_from"])

        self.interface_input_from = None
        self.interface_input_to = None
        self.interface_output_to = None

        # run time
        self.run_time = 0.0

    def initialize(self, interface_input_from, interface_output_to):
        super().initialize()

        self.solver_wrapper.initialize()

        # create input mapper
        self.interface_input_from = interface_input_from.copy()
        self.interface_input_to = self.solver_wrapper.get_interface_input()
        self.mapper_interface_input.initialize(self.interface_input_from, self.interface_input_to)

        # create output mapper
        self.interface_output_to = interface_output_to.copy()
        self.interface_output_from = self.solver_wrapper.get_interface_output()
        self.mapper_interface_output.initialize(self.interface_output_from, self.interface_output_to)

        #create mp input_update_to
        model1 = self.interface_input_to.model
        parameters1 = self.interface_input_to.parameters
        for dct in parameters1:
            dct['variables'] = ["displacement"]
        self.new_interface_input_to = data_structure.Interface(parameters1,model1)

        #create mp_input_udate_from
        model2 = self.interface_input_from.model
        parameters2 = self.interface_input_from.parameters
        for dct2 in parameters2:
            dct2['variables'] = ["displacement"]
        self.new_interface_input_from = data_structure.Interface(parameters2, model2)

        for i in self.interface_output_to.parameters:
            tmp = self.interface_output_to.get_model_part(i['model_part'])
            print(tmp)
            # for j in range(tmp.x0.size):

                # print(tmp.x0[j])
                # print(tmp.y0[j])
                # print(tmp.z0[j])
            self.coords_in = np.column_stack((tmp.x0, tmp.y0, tmp.z0))
            print(self.coords_in)
        #create input_update_to mapper
        self.mapper_interface_input_update_to.initialize(self.interface_output_from, self.new_interface_input_to)


        #create input_update_from mapper
        self.mapper_interface_input_update_from.initialize(self.interface_output_to, self.new_interface_input_from)

    def initialize_solution_step(self):
        super().initialize_solution_step()

        self.solver_wrapper.initialize_solution_step()
        self.iteration = 0

    @tools.time_solve_solution_step
    def solve_solution_step(self, interface_input_from):

        self.iteration += 1

        self.interface_input_from = interface_input_from.copy()

        # Creating an updated self.interface_input_from

        # purple box 1 mapping
        self.mapper_interface_input_update_from(self.interface_output_to, self.new_interface_input_from)

        # add the mapped displacement to interface_input_from by creating a new mp
        self.model_new_from = data_structure.Model()
        for item_new_input_from in self.new_interface_input_from.parameters:
            mp_intermediate_input_from = self.new_interface_input_from.get_model_part(item_new_input_from['model_part'])
            mp_intermediate_from = self.new_interface_input_from.get_variable_data(item_new_input_from['model_part'],
                                                                               item_new_input_from['variables'][0])
            x0_new = mp_intermediate_from[:, 0] + mp_intermediate_input_from.x0
            y0_new = mp_intermediate_from[:, 1] + mp_intermediate_input_from.y0
            z0_new = mp_intermediate_from[:, 2] + mp_intermediate_input_from.z0
            ids_from = np.arange(0, mp_intermediate_input_from.x0.size)

            # create new mp with interface_new_to to update self.interface_input_to
            self.model_new_from.create_model_part(item_new_input_from['model_part'],
                                                x0_new, y0_new, z0_new,
                                                ids_from)
            self.parameters_new_from = [{'model_part': item_new_input_from['model_part'],
                                  'variables': ['pressure', 'traction']}]

        self.interface_new_from = data_structure.Interface(self.parameters_new_from, self.model_new_from)

        tmp = {}
        for item_input_from in self.interface_input_from.parameters:
            j = 0
            for i in item_input_from['variables']:
                self.tmp_[j] = self.interface_input_from.get_variable_data(item_input_from['model_part'], item_input_from['variables'][j])
                self.interface_new_from.set_variable_data(item_new_input_from['model_part'],
                                                    ['variables'][j], self.tmp_[j])
                j+=1


        self.interface_input_from = self.interface_new_from

        # Creating an updated self.interface_input_to

        #purple box 2 mapping
        self.mapper_interface_input_update_to(self.interface_output_from, self.new_interface_input_to)

        # add the mapped displacement to interface_input_to by creating a new mp
        self.model_new_to = data_structure.Model()
        for item_new_input_to in self.new_interface_input_to.parameters:
            mp_intermediate_input_to = self.new_interface_input_to(item_new_input_to['model_part'])
            tmp = self.new_interface_input_to.get_variable_data(item_new_input_to['model_part'],
                                                                               item_new_input_to['variables'][0])
            x0_new_to = tmp[:, 0] + mp_intermediate_input_to.x0
            y0_new_to = tmp[:, 1] + mp_intermediate_input_to.y0
            z0_new_to = tmp[:, 2] + mp_intermediate_input_to.z0
            ids_to = np.arange(0, mp_intermediate_input_to.x0.size)

            #create new mp with interface_new_to to update self.interface_input_to
            self.model_new_to.create_model_part(item_new_input_to['model_part'],
                                                x0_new_to, y0_new_to, z0_new_to,
                                                ids_to)
            self.parameters_new_to = [{'model_part': item_new_input_to['model_part'],
                                  'variables': ['pressure', 'traction']}]

        self.interface_new_to = data_structure.Interface(self.parameters_new_to, self.model_new_to)

        self.interface_input_to = self.interface_new_to

        self.mapper_interface_input(self.interface_input_from, self.interface_input_to)

        interface_output_from = self.solver_wrapper.solve_solution_step(self.interface_input_to)

        # Creating an updated interface_output_from

        self.model_new_output_from = data_structure.Model()
        #get the point displacements of current calculation in the structural solver
        for item_output_from in self.interface_output_from.parameters:
            mp_output_from = self.interface_output_from.get_model_part(item_output_from['model_part'])
            varia = self.interface_output_from.get_variable_data(item_output_from['model_part'], item_output_from['variables'][0])
            x0_new_from = np.add(varia[:,0],mp_output_from.x0)
            y0_new_from = np.add(varia[:,1],mp_output_from.y0)
            z0_new_from = np.add(varia[:,2],mp_output_from.z0)
            ids = np.arange(0, mp_output_from.x0.size)

            #creating new mp
            self.model_new_output_from.create_model_part( item_output_from['model_part'], x0_new_from, y0_new_from, z0_new_from, ids)
            self.parameters_new_output_from = [{'model_part': item_output_from['model_part'], 'variables' :  item_output_from['variables'][0] }]

        self.interface_output_from_new = data_structure.Interface(self.parameters_new_output_from, self.model_new_output_from)
        self.interface_output_from_new.set_variable_data(item_output_from['model_part']+str(self.iteration),item_output_from['variables'][0], varia)
        self.interface_output_from = self.interface_output_from_new.get_model_part(item_output_from['model_part']+str(self.iteration))

         #store the displacements for following iteration

        self.interface_output_from = self.interface_output_from.copy()

        self.mapper_interface_output(interface_output_from, self.interface_output_to)

        return self.interface_output_to.copy()

    def finalize_solution_step(self):
        super().finalize_solution_step()

        self.solver_wrapper.finalize_solution_step()

    def finalize(self):
        super().finalize()

        self.solver_wrapper.finalize()
        self.mapper_interface_input.finalize()
        self.mapper_interface_output.finalize()
        self.mapper_interface_input_update_from.finalize()
        self.mapper_interface_input_update_to.finalize()

    def output_solution_step(self):
        super().output_solution_step()

        self.solver_wrapper.output_solution_step()
        self.mapper_interface_input.output_solution_step()
        self.mapper_interface_output.output_solution_step()
        self.mapper_interface_input_update_from.output_solution_step()
        self.mapper_interface_input_update_to.output_solution_step()

    def get_interface_input(self):
        # does not contain most recent data
        return self.interface_input_from

    def get_interface_output(self):
        interface_output_from = self.solver_wrapper.get_interface_output()
        self.mapper_interface_output(interface_output_from, self.interface_output_to)
        return self.interface_output_to

    def print_components_info(self, pre):
        tools.print_info(pre, "The component ", self.__class__.__name__, " maps the following solver wrapper:")
        pre = tools.update_pre(pre)
        self.solver_wrapper.print_components_info(pre + '├─')
        tools.print_info(pre, '├─', "Input mapper:")
        self.mapper_interface_input.print_components_info(pre + '│ └─')
        tools.print_info(pre, '└─', "Output mapper:")
        self.mapper_interface_output.print_components_info(pre + '  └─')
