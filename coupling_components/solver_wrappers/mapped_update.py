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
        self.model_to = data_structure.Model()

        for item_input_to in self.interface_input_to.parameters:
            mp_input_to = self.interface_input_to.get_model_part(item_input_to['model_part'])
            self.model_to.create_model_part('new_mp_input_to', mp_input_to.x0, mp_input_to.y0, mp_input_to.z0, np.arange(mp_input_to.x0.size))
            parameters_to = [{'model_part':'new_mp_input_to', 'variables':['displacement']}]
            self.new_interface_input_to = data_structure.Interface(parameters_to, self.model_to)
            self.new_interface_input_to.set_variable_data('new_mp_input_to', 'displacement', np.zeros(mp_input_to.x0.size,3))

        #create mp_input_udate_from
        self.model_from = data_structure.Model()

        for item_input_from in self.interface_input_from.parameters:
            mp_input_from = self.interface_input_from.get_model_part(item_input_from['model_part'])
            self.model_from.create_model_part('new_mp_input_from', mp_input_from.x0, mp_input_from.y0, mp_input_from.z0, np.arange(mp_input_from.x0.size))
            parameters_from = [{'model_part':'new_mp_input_from', 'variables':['displacement']}]
            self.new_interface_input_from = data_structure.Interface(parameters_from, self.model_from)
            self.new_interface_input_from.set_variable_data('new_mp_input_from', 'displacement',np.zeros(mp_input_from.x0.size, 3))

        #create input_update_to mapper
        self.mapper_interface_input_update_to.initialize(self.interface_output_from, self.new_interface_input_to)

        #create input_update_from mapper
        self.mapper_interface_input_update_from.initialize(self.interface_output_to, self.new_interface_input_from)

    def initialize_solution_step(self):
        super().initialize_solution_step()

        self.solver_wrapper.initialize_solution_step()
        self.iteration = 0
        self.update_array_to =np.zeros((self.interface_output_from.x0.size,1))
        self.update_array_from = np.zeros((self.interface_output_to.x0.size, 1))

    @tools.time_solve_solution_step
    def solve_solution_step(self, interface_input_from):

        self.iteration += 1

        self.interface_input_from = interface_input_from.copy()

        # Creating an updated self.interface_input_from

        # create an interface from output_to with the displacements of previous iteration
        self.interface_output_to_intermediate = self.interface_output_to.copy()
        for item_output_to in self.interface_output_to.parameters:
            self.interface_output_to_intermediate.set_variable_data(item_output_to['model_part'], item_output_to
                ['variables'][0], self.update_array_from)

        # purple box 1 mapping
        self.mapper_interface_input_update_from(self.interface_output_to_intermediate, self.new_interface_input_from)

        # add the mapped displacement to interface_input_from by creating a new mp
        self.model_new_from = data_structure.Model()
        for item_new_input_from in self.new_interface_input_from.parameters:
            mp_intermediate_input_from = self.new_interface_input_from(item_new_input_from['model_part'])
            mp_intermediate_from = self.new_interface_input_from.get_variable_data(item_new_input_from['model_part'],
                                                                               item_new_input_from['variables'][0])
            new_mp_input_from = np.zeros((mp_intermediate_input_from.x0.size, 3))
            new_mp_input_from[:, 0] = np.add(mp_intermediate_from[:, 0], mp_intermediate_input_from.x0)
            new_mp_input_from[:, 1] = np.add(mp_intermediate_from[:, 1], mp_intermediate_input_from.y0)
            new_mp_input_from[:, 2] = np.add(mp_intermediate_from[:, 2], mp_intermediate_input_from.z0)
            ids_from = np.arange(0, mp_intermediate_input_from.x0.size)

            # create new mp with interface_new_to to update self.interface_input_to
            self.model_new_to.create_model_part(item_new_input_from['model_part'] + str(self.iteration),
                                                new_mp_input_from[:, 0], new_mp_input_from[:, 1], new_mp_input_from[:, 2],
                                                ids_from)
            parameters_new_from = [{'model_part': item_new_input_from['model_part'] + str(self.iteration),
                                  'variables': ['pressure', 'traction']}]
            self.interface_new_from = data_structure.Interface(parameters_new_from, self.model_new_from)
            for item_input_from in self.interface_input_from.parameters:
                tmp_pressure = self.interface_input_from.get_variable_data(item_input_from['model_part'], item_input_from['variables'][0])
                tmp_traction = self.interface_input_from.get_variable_data(item_input_from['model_part'], item_input_from['variables'][1])

            self.interface_new_from.set_variable_data(item_new_input_from['model_part'] + str(self.iteration),
                                                    ['variables'][0], tmp_pressure)
            self.interface_new_from.set_variable_data(item_new_input_from['model_part'] + str(self.iteration),
                                                    ['variables'][1], tmp_traction)

        self.interface_input_from = self.interface_new_from
        # Creating an updated self.interface_input_to

        # create an interface from output_from with the displacements of previous iteration
        self.interface_output_from_intermediate = self.interface_output_from.copy()
        for item_output_from in self.interface_output_from.parameters:
            self.interface_output_from_intermediate.set_variable_data(item_output_from['model_part'], item_output_from
                ['variables'][0], self.update_array_to)

        #purple box 2 mapping
        self.mapper_interface_input_update_to(self.interface_output_from_intermediate, self.new_interface_input_to)

        # add the mapped displacement to interface_input_to by creating a new mp
        self.model_new_to = data_structure.Model()
        for item_new_input_to in self.new_interface_input_to.parameters:
            mp_intermediate_input_to = self.new_interface_input_to(item_new_input_to['model_part'])
            mp_intermediate_to = self.new_interface_input_to.get_variable_data(item_new_input_to['model_part'],
                                                                               item_new_input_to['variables'][0])
            new_mp_input_to = np.zeros((mp_intermediate_input_to.x0.size, 3))
            new_mp_input_to[:, 0] = np.add(mp_intermediate_to[:, 0], mp_intermediate_input_to.x0)
            new_mp_input_to[:, 1] = np.add(mp_intermediate_to[:, 1], mp_intermediate_input_to.y0)
            new_mp_input_to[:, 2] = np.add(mp_intermediate_to[:, 2], mp_intermediate_input_to.z0)
            ids_to = np.arange(0, mp_intermediate_input_to.x0.size)

            #create new mp with interface_new_to to update self.interface_input_to
            self.model_new_to.create_model_part(item_new_input_to['model_part'] + str(self.iteration),
                                                new_mp_input_to[:, 0], new_mp_input_to[:, 1], new_mp_input_to[:, 2],
                                                ids_to)
            parameters_new_to = [{'model_part': item_new_input_to['model_part'] + str(self.iteration),
                                  'variables': ['pressure', 'traction']}]
            self.interface_new_to = data_structure.Interface(parameters_new_to, self.model_new_to)
            self.interface_new_to.set_variable_data(item_new_input_to['model_part'] + str(self.iteration),
                                                    ['variables'][0], np.zeros((mp_intermediate_input_to.x0.size,1)))
            self.interface_new_to.set_variable_data(item_new_input_to['model_part'] + str(self.iteration),
                                                    ['variables'][1], np.zeros((mp_intermediate_input_to.x0.size, 3)))

        self.interface_input_to = self.interface_new_to

        self.mapper_interface_input(self.interface_input_from, self.interface_input_to)

        interface_output_from = self.solver_wrapper.solve_solution_step(self.interface_input_to)

        # Creating an updated interface_output_from

        self.model_new_output_from = data_structure.Model()
        #get the point displacements of current calculation in the structural solver
        for item_output_from in self.interface_output_from.parameters:
            mp_output_from = self.interface_output_from.get_model_part(item_output_from['model_part'])
            varia = self.interface_output_from.get_variable_data(item_output_from['model_part'], item_output_from['variables'][0])
            new_mp = np.zeros((mp_output_from.x0.size,3))
            new_mp[:,0] = np.add(varia[:,0],mp_output_from.x0)
            new_mp[:,1] = np.add(varia[:,1],mp_output_from.y0)
            new_mp[:,2] = np.add(varia[:,2],mp_output_from.z0)
            ids = np.arange(0, mp_output_from.x0.size)

            #creating new mp
            self.model_new_output_from.create_model_part( item_output_from['model_part']+str(self.iteration), new_mp[:, 0], new_mp[:, 1], new_mp[:, 2], ids)
            parameters = [{'model_part': item_output_from['model_part']+str(self.iteration), 'variables' :  item_output_from['variables'][0] }]
            self.interface_output_from_new = data_structure.Interface(parameters, self.model_new_output_from)
            self.interface_output_from_new.set_variable_data(item_output_from['model_part']+str(self.iteration),item_output_from['variables'][0], varia)
            self.interface_output_from = self.interface_output_from_new.get_model_part(item_output_from['model_part']+str(self.iteration))

            #store the displacements for following iteration
            self.update_array_to = varia

        self.mapper_interface_output(interface_output_from, self.interface_output_to)

        # store the displacements for following iteration
        for item_output_to in self.interface_output_to.parameters:
            varia2 = self.interface_output_to.get_variable_data(item_output_to['model_part'], item_output_to['variables'][0])
            self.update_array_from = varia2

        return self.interface_output_to

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
