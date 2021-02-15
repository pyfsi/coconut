from coconut.coupling_components.tools import create_instance
from coconut.coupling_components.component import Component
from coconut.coupling_components import tools
import numpy as np
from coconut import data_structure

""" proposed changes to mapped.py
- do initialization of mappers in Initialize method, would be more logical
- remove all set_interface_input/Output methods?
- use copy in get_interface_input/Output methods?
    and just refer to actual solver wrapper in SolverWrapperMapped
- all Interfaces are stored in this mapper, e.g. self.interface_output_to and 3 others;
    I see no reason for this; furthermore, it is only useful to store it if you take copies all the time
- output_solution_step is barely used; what's the deal with it??
"""


def create(parameters):
    return SolverWrapperMapped(parameters)


class SolverWrapperMapped(Component):
    def __init__(self, parameters):
        super().__init__()

        # Read parameters
        self.parameters = parameters
        self.settings = parameters["settings"]

        # Create solver
        self.solver_wrapper = create_instance(self.settings["solver_wrapper"])

        # run time
        self.run_time = 0.0

    def initialize(self):
        super().initialize()

        self.solver_wrapper.initialize()

    def initialize_solution_step(self):
        super().initialize_solution_step()

        self.solver_wrapper.initialize_solution_step()
        self.iteration = 0
        self.updateArray = {}

    @tools.time_solve_solution_step
    def solve_solution_step(self, interface_input_from):
        self.iteration += 1
        self.interface_input_from = interface_input_from

        for item_input_to in self.interface_input_to.parameters:
            #create MP with displacement field in input_to mapping
            mp_input_to = self.interface_input_to.get_model_part(item_input_to['model_part'])
            self.model.create_model_part('intermediate_mp' + str(self.iteration), mp_input_to.x0, mp_input_to.y0,
                                         mp_input_to.z0,
                                         np.arange(mp_input_to.x0.size))
            parameters = [{'model_part': 'intermediate_mp' + str(self.iteration), 'variables': ['displacement']}]
            self.interface_intermediate = data_structure.Interface(parameters, self.model)
            self.interface_intermediate.set_variable_data('intermediate_mp' + str(self.iteration), 'displacement',
                                                          np.zeros((mp_input_to.x0.size, 3)))
            mp_intermediate = self.interface_intermediate.get_variable_data('intermediate_mp' + str(self.iteration),
                                                                            'displacement')
            print("zeros")
            print(mp_intermediate)
            #create a variable with the displcaments of previous iteration
            if self.iteration > 1:
                self.interface_output_from_intermediate = self.interface_output_from.copy()
                self.interface_output_from_intermediate.set_variable_data('BEAMINSIDEMOVING_nodes', 'displacement',
                                                                          self.myArrays[self.iteration - 1])
                check = self.interface_output_from_intermediate.get_variable_data('BEAMINSIDEMOVING_nodes',
                                                                                  'displacement')
                self.mapper_interface_inputUpdate(self.interface_output_from_intermediate, self.interface_intermediate)
                # add the mapped displacement to interface_input_to by creating new MP
                for item_input_intermediat_to in self.interface_intermediate:
                    mp_intermediate_to = self.interface_intermediate(item_input_intermediat_to['model_part'])
                    mp_intermediate = self.interface_intermediate.get_variable_data('intermediate_mp' + str(self.iteration),
                                                                                'displacement')
                    new_mp_input_to = [mp_intermediate_to.x0.size, 3]
                    new_mp_input_to[:,0] = np.add(mp_intermediate[:,0],mp_intermediate_to.x0)
                    new_mp_input_to[:,1] = np.add(mp_intermediate[:,1],mp_intermediate_to.y0)
                    ids_to = np.arange(0,mp_intermediate_to.x0.size)
                    self.model.create_model_part('new_to' +str(self.iteration), new_mp_input_to[:,0], new_mp_input_to[:,1], new_mp_input_to[:,2], ids_to)
                    self.interface_input_to = self.model.get_model_part('new_to' + str(self.iteration))


        self.mapper_interface_input(self.interface_input_from, self.interface_input_to)

        self.interface_output_from = self.solver_wrapper.solve_solution_step(self.interface_input_to)
        for item_output_from in (self.interface_output_from.parameters):
            mp_output_from = self.interface_output_from.get_model_part(item_output_from['model_part'])
            varia = self.interface_output_from.get_variable_data(item_output_from['model_part'], item_output_from['variables'][0])
            new_mp = [mp_output_from.X0.size, 3]
            new_mp[:,0] = np.add(varia[:,0], mp_output_from.x0)
            new_mp[:,1] = np.add(varia[:,1], mp_output_from.y0)
            ids = np.arange(0, mp_output_from.X0.size)
            self.model.create_model_part('new_mp', new_mp[:, 0], new_mp[:, 1], new_mp[:, 2], ids)
            self.interface_output_from = self.model.get_model_part('new_mp')
            # save the displacement for transfer towards self.interface_input_to
            self.updateArray[self.iteration] = varia

        self.mapper_interface_output(self.interface_output_from, self.interface_output_to)

        return self.interface_output_to

    def finalize_solution_step(self):
        super().finalize_solution_step()

        self.solver_wrapper.finalize_solution_step()

    def finalize(self):
        super().finalize()

        self.solver_wrapper.finalize()
        self.mapper_interface_input.finalize()
        self.mapper_interface_output.finalize()

    def output_solution_step(self):
        super().output_solution_step()

        self.solver_wrapper.output_solution_step()
        self.mapper_interface_input.output_solution_step()
        self.mapper_interface_output.output_solution_step()

    def get_interface_input(self):
        # Does not contain most recent data
        # *** shouldn't this just call the underlying solver wrapper?
        return self.interface_input_from

    def set_interface_input(self, interface_input_from):
        # Create input mapper
        self.interface_input_from = interface_input_from.copy()
        self.interface_input_to = self.solver_wrapper.get_interface_input()

        self.mapper_interface_input = create_instance(self.settings["mapper_interface_input"])
        self.mapper_interface_input.initialize(self.interface_input_from, self.interface_input_to)

    def set_interface_inputUpdate_to(self, interface_output_from):
        # Create input mapper
        self.interface_output_from = interface_output_from.copy()
        self.interface_input_to = self.solver_wrapper.get_interface_output()

        self.mapper_interface_input = create_instance(self.settings["mapper_interface_inputUpdate_to"])
        self.mapper_interface_input.initialize(self.interface_output_from, self.interface_input_to)

    def get_interface_output(self):
        self.interface_output_from = self.solver_wrapper.get_interface_output()
        self.mapper_interface_output(self.interface_output_from, self.interface_output_to)
        return self.interface_output_to.copy()

    def set_interface_output(self, interface_output_to):
        # Create output mapper
        self.interface_output_to = interface_output_to.copy()
        self.interface_output_from = self.solver_wrapper.get_interface_output()

        self.mapper_interface_output = create_instance(self.settings["mapper_interface_output"])
        self.mapper_interface_output.initialize(self.interface_output_from, self.interface_output_to)

    def print_components_info(self, pre):
        tools.print_info(pre, "The component ", self.__class__.__name__, " maps the following solver wrapper:")
        pre = tools.update_pre(pre)
        self.solver_wrapper.print_components_info(pre + '├─')
        tools.print_info(pre, '├─', "Input mapper:")
        self.mapper_interface_input.print_components_info(pre + '│ └─')
        tools.print_info(pre, '└─', "Output mapper:")
        self.mapper_interface_output.print_components_info(pre + '  └─')
