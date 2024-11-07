from coconut.coupling_components.solver_wrappers.solver_wrapper import SolverWrapper
from coconut import tools
from coconut.data_structure import Model, Interface
import coconut.coupling_components.solver_wrappers.python.banded as bnd

import numpy as np
import os
from os.path import join
from scipy.linalg import solve_banded
import json
import pickle


def create(parameters):
    return SolverWrapperSaturatedSolid(parameters)


class SolverWrapperSaturatedSolid(SolverWrapper):
    check_coupling_convergence_possible = False  # can solver check convergence after 1 iteration?

    @tools.time_initialize
    def __init__(self, parameters):
        super().__init__(parameters)

        # reading
        self.parameters = parameters
        self.settings = parameters['settings']
        self.working_directory = self.settings['working_directory']
        input_file = self.settings.get('input_file')
        if input_file is not None:
            case_file_name = join(self.working_directory, input_file)
            with open(case_file_name, 'r') as case_file:
                case_file_settings = json.load(case_file)
            case_file_settings.update(self.settings)
            with open(case_file_name, 'w') as case_file:
                json.dump(case_file_settings, case_file, indent=2)
            self.settings.update(case_file_settings)

        # store settings
        self.dt = self.settings['delta_t']
        self.restart = self.settings['save_restart']
        self.timestep = self.settings['timestep_start']
        self.interface_settings = self.settings['interface']
        self.material_properties = self.settings['material_properties']
        self.mapper_settings = self.settings['conservative_mapper']

        if self.restart != 0:  # restart
            raise NotImplementedError('Restart functionality not supported.')

        self.x0 = self.interface_settings['x0']
        self.y0 = self.interface_settings['y0']
        self.x1 = self.interface_settings['x1']
        self.y1 = self.interface_settings['y1']
        self.n = self.interface_settings['faces']
        self.mov_dir = self.interface_settings['movement_direction']
        self.mov_dir = np.array((self.mov_dir))
        
        self.rho = self.material_properties['rho']
        self.L = self.material_properties['latent']

        self.mapping_domain = self.mapper_settings['mapping_domain']
        self.mapping_limits = self.mapper_settings['mapping_limits']
        self.conservative = True

        # initialization
        x = np.linspace(self.x0, self.x1, self.n + 1)
        y = np.linspace(self.y0, self.y1, self.n + 1)
        ini_coord_faces = []
        for i in range(self.n):
            ini_coord_faces.append([(x[i] + x[i+1]) / 2, (y[i] + y[i+1]) / 2, 0])
        self.ini_coord_faces = np.array(ini_coord_faces)

        x = np.reshape(x, (self.n + 1, 1))
        y = np.reshape(y, (self.n + 1, 1))
        z = np.zeros((self.n + 1, 1))
        self.ini_coord_nodes = np.concatenate((x, y), axis=1)
        self.ini_coord_nodes = np.concatenate((self.ini_coord_nodes, z), axis=1)
        self.area = self.area_calc(self.ini_coord_nodes)

        self.prev_face_disp = np.zeros((self.n, 3))  # previous total face displacement
        self.prev_disp = np.zeros((self.n + 1, 3))  # previous total node displacement

        self.face_dx = np.zeros((self.n, 3))  # latest time step face displacement
        self.dx = np.zeros((self.n + 1, 3))  # latest time step node displacement

        self.heat_flux = np.zeros((self.n,1)) # heat flux [W/m^2]

        # create input & output ModelParts
        self.model = Model()

        self.input_mp_name = self.settings['interface_input'][0]['model_part']
        self.output_mp_name = self.settings['interface_output'][0]['model_part']

        self.model.create_model_part(self.input_mp_name, self.ini_coord_faces[:,0].flatten(),
                                               self.ini_coord_faces[:,1].flatten(), np.zeros(self.n), np.arange(self.n))
        self.model.create_model_part(self.output_mp_name, self.ini_coord_nodes[:,0].flatten(),
                                               self.ini_coord_nodes[:,1].flatten(), np.zeros(self.n + 1), np.arange(self.n + 1))

        # input & output interfaces
        self.interface_input = Interface(self.settings['interface_input'], self.model)
        self.interface_input.set_variable_data(self.input_mp_name, 'heat_flux', self.heat_flux)

        self.interface_output = Interface(self.settings['interface_output'], self.model)
        self.interface_output.set_variable_data(self.output_mp_name, 'displacement', self.prev_disp + self.dx)

        # internal models for face-to-node mapping
        self.internal_model = Model()
        self.internal_face_settings = [{"model_part": self.output_mp_name, "variables": ["1ts_disp", "prev_disp", "area"]}]
        self.internal_node_settings = [{"model_part": self.output_mp_name, "variables": ["displacement", "1ts_disp", "prev_disp"]}]

        self.internal_model.create_model_part(self.output_mp_name, self.ini_coord_faces[:, 0].flatten(),
                                     self.ini_coord_faces[:, 1].flatten(), np.zeros(self.n), np.arange(self.n))

        self.interface_internal_faces = Interface(self.internal_face_settings, self.internal_model)
        self.interface_internal_nodes = Interface(self.internal_node_settings, self.model)

        # create and initialize face to node (f2n) displacement mapper
        f2n_settings = {"type": "mappers.interface", "settings": {"type": "mappers.linear_conservative",
                                                                  "settings": {"directions": ["x", "y"],
                                                                               "check_bounding_box": False,
                                                                               "domain": self.mapping_domain,
                                                                               "limits": self.mapping_limits}}}
        self.mapper_f2n = tools.create_instance(f2n_settings)
        self.mapper_f2n.initialize(self.interface_internal_faces, self.interface_internal_nodes)

        # create and initialize node to face (n2f) displacement mapper
        n2f_settings = {"type": "mappers.interface", "settings": {"type": "mappers.linear",
                                                                  "settings": {"directions": ["x", "y"],
                                                                               "check_bounding_box": False}}}
        self.mapper_n2f = tools.create_instance(n2f_settings)
        self.mapper_n2f.initialize(self.interface_internal_nodes, self.interface_internal_faces)

        # create and initialize node to node (n2n) displacement mapper for output purposes
        n2n_settings = {"type": "mappers.interface", "settings": {"type": "mappers.nearest",
                                                                  "settings": {"directions": ["x", "y"],
                                                                               "check_bounding_box": False}}}
        self.mapper_n2n = tools.create_instance(n2n_settings)
        self.mapper_n2n.initialize(self.interface_internal_nodes, self.interface_output)

    @tools.time_initialize
    def initialize(self):
        super().initialize()

    def initialize_solution_step(self):
        super().initialize_solution_step()

        self.timestep += 1
        self.mapper_n2f.map_n2f(self.interface_internal_nodes, self.interface_internal_faces, conservative=self.conservative)

        self.prev_disp += self.dx
        self.interface_internal_nodes.set_variable_data(self.output_mp_name, 'prev_disp', self.prev_disp)

    @tools.time_solve_solution_step
    def solve_solution_step(self, interface_input):
        # input
        self.interface_input = interface_input.copy()
        self.heat_flux = interface_input.get_variable_data(self.input_mp_name, 'heat_flux') # [W/m^2]
        self.heat_flux = -1*self.heat_flux # negative heat flux for liquid domain is positive for solid domain
        
        disp_magn = (self.heat_flux * self.dt) / (self.rho * self.L) # Stefan condition
        self.area, normal_array = self.area_calc(self.ini_coord_nodes + self.prev_disp)

        self.face_dx = disp_magn*normal_array

        self.interface_internal_faces.set_variable_data(self.output_mp_name, '1ts_disp', self.face_dx)
        self.interface_internal_faces.set_variable_data(self.output_mp_name, 'area', self.area)

        # map face displacement to node displacement
        self.mapper_f2n.map_f2n(self.interface_internal_faces, self.interface_internal_nodes, conservative=self.conservative)
        
        self.dx = self.interface_internal_nodes.get_variable_data(self.output_mp_name, '1ts_disp')
        self.interface_output.set_variable_data(self.output_mp_name, 'displacement', self.prev_disp + self.dx)
        self.interface_internal_nodes.set_variable_data(self.output_mp_name, 'displacement', self.prev_disp + self.dx)

        # output
        return self.interface_output

    def finalize_solution_step(self):
        super().finalize_solution_step()

    def finalize(self):
        super().finalize()
    
    def area_calc(self, nodes):
        """
        Receives: array with ordered (n+1) node coordinates
        Returns: array with ordered (n) face areas bounded by given nodes
        """

        nf = np.shape(nodes)[0] - 1
        area = np.zeros(nf)
        normal_array = np.zeros((nf, 3))

        for i in range(nf):
            area[i] = np.sqrt((nodes[i+1][0] - nodes[i][0])**2 + (nodes[i+1][1] - nodes[i][1])**2)
            nx = (nodes[i+1][0] - nodes[i][0]) / area[i]
            ny = (nodes[i + 1][1] - nodes[i][1]) / area[i]
            normal = np.array([ny, -1*nx])
            if np.dot(self.mov_dir, normal) < 0:
                normal = -1*normal
            normal_array[i, 0:2] = normal

        return np.reshape(area, (nf, 1)), normal_array