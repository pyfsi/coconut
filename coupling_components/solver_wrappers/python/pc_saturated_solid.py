from coconut.coupling_components.solver_wrappers.solver_wrapper import SolverWrapper
from coconut import tools
from coconut.data_structure import Model, Interface
import coconut.coupling_components.solver_wrappers.python.banded as bnd

import numpy as np
import os
from os.path import join
from scipy.linalg import solve_banded
from scipy.spatial import cKDTree
import json
import pickle
import pandas as pd


def create(parameters):
    return SolverWrapperSaturatedSolid(parameters)


class SolverWrapperSaturatedSolid(SolverWrapper):
    check_coupling_convergence_possible = False  # can solver check convergence after 1 iteration?

    # define input and output variables
    accepted_in_var = ['displacement', 'heat_flux']
    accepted_out_var = ['displacement']

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
        self.dt = self.settings['timestep_size']
        self.timestep_start = self.settings.get('timestep_start', 0)
        self.timestep = self.timestep_start  # time step
        self.save_restart = self.settings.get('save_restart', 0)
        self.interface_settings = self.settings['interface']
        self.material_properties = self.settings['material_properties']
        self.mapper_settings = self.settings['conservative_mapper']

        self.mov_dir = self.interface_settings['movement_direction']
        # Check for string-based direction ('inward'/'outward') or legacy vector
        if isinstance(self.mov_dir, str):
            if self.mov_dir not in ['inward', 'outward']:
                raise ValueError(f"Invalid 'movement_direction': '{self.mov_dir}'. Must be 'inward' or 'outward'.")
        else:
            self.mov_dir = np.array(self.mov_dir)
        self.closed = self.interface_settings.get('closed', False)

        self.rho = self.material_properties['rho']
        self.L = self.material_properties['latent']

        self.mapping_domain = self.mapper_settings['mapping_domain']
        self.mapping_limits = self.mapper_settings['mapping_limits']
        self.conservative = self.mapper_settings.get('mapping_conservative', False)
        self.vol_tolerance = self.mapper_settings.get('volume_tolerance', 1e-14)
        self.mapper_iterations = self.mapper_settings.get('mapper_iterations', 20)
        self.mapper_relaxation = self.mapper_settings.get('mapper_relaxation', 0.4)
        self.mapper_projection = self.mapper_settings.get('projection_order', 1) # 0: no projection, 1: 1st order upwind projection, 2: 2nd order upwind

        # initialization
        interface_file = self.interface_settings.get('interface_file')
        if self.timestep_start == 0:  # no restart
            if interface_file is not None:
                interface_path = join(self.working_directory, interface_file)
                try:
                    df_interface = pd.read_csv(interface_path)
                except FileNotFoundError:
                    raise FileNotFoundError(f"The specified interface file was not found: {interface_path}")

                raw_x = df_interface['x-coordinate'].values
                raw_y = df_interface['y-coordinate'].values
                points = np.column_stack((raw_x, raw_y))

                # --- OPTIMIZED SORTING (KD-Tree) ---
                if len(points) > 1:
                    # 1. Build the KD-Tree (Very fast spatial index)
                    tree = cKDTree(points)

                    # 2. Pick a starting point.
                    current_idx = np.argmin(points[:, 0])

                    sorted_indices = [current_idx]
                    visited = set([current_idx])

                    # 3. Traverse the chain
                    for _ in range(len(points) - 1):
                        # Query k=2 nearest neighbors (1st is the point itself, 2nd is the neighbor)
                        # We query k=5 just in case the immediate neighbors are already visited
                        dists, indices = tree.query(points[current_idx], k=min(10, len(points)))

                        found_next = False
                        for neighbor_idx in indices:
                            if neighbor_idx not in visited:
                                visited.add(neighbor_idx)
                                sorted_indices.append(neighbor_idx)
                                current_idx = neighbor_idx
                                found_next = True
                                break

                        if not found_next:
                            # This happens if the surface is disjoint or you hit a dead end
                            # but unvisited points remain elsewhere.
                            break

                    # Apply the sort
                    points = points[sorted_indices]

                x_nodes = points[:, 0]
                y_nodes = points[:, 1]

                self.nn = len(x_nodes)
                self.nf = self.nn if self.closed else self.nn - 1

                print(f'Solid nodes: {self.nn}')
                print(f'Solid faces: {self.nf}')

                z_nodes = np.zeros(self.nn)
                self.ini_coord_nodes = np.vstack((x_nodes, y_nodes, z_nodes)).T

                # Calculate initial face coordinates as midpoints of the nodes
                self.ini_coord_faces = (self.ini_coord_nodes + np.roll(self.ini_coord_nodes, -1, axis=0)) / 2
                if not self.closed:
                    self.ini_coord_faces = self.ini_coord_faces[:-1]
            else:  # Fallback to original straight line logic if file is not provided
                if not self.closed:
                    self.x0 = self.interface_settings['x0']
                    self.y0 = self.interface_settings['y0']
                    self.x1 = self.interface_settings['x1']
                    self.y1 = self.interface_settings['y1']
                    self.nf = self.interface_settings['faces']
    
                    x = np.linspace(self.x0, self.x1, self.nn)
                    y = np.linspace(self.y0, self.y1, self.nn)
                    ini_coord_faces = []
                    for i in range(self.nf):
                        ini_coord_faces.append([(x[i] + x[i+1]) / 2, (y[i] + y[i+1]) / 2, 0])
                    self.ini_coord_faces = np.array(ini_coord_faces) # initial face coordinates
    
                    x = np.reshape(x, (self.nn, 1))
                    y = np.reshape(y, (self.nn, 1))
                    z = np.zeros((self.nn, 1))
                    self.ini_coord_nodes = np.concatenate((x, y, z), axis=1)

            self.prev_face_disp = np.zeros((self.nf, 3))  # previous total face displacement
            self.prev_disp = np.zeros((self.nn, 3))  # previous total node displacement

            self.face_dx = np.zeros((self.nf, 3))  # latest time step face displacement
            self.dx = np.zeros((self.nn, 3))  # latest time step node displacement
        else: # restart
            file_name = join(self.working_directory, f'case_timestep{self.timestep_start}.pickle')
            with open(file_name, 'rb') as file:
                data = pickle.load(file)
            self.ini_coord_nodes = data['ini_nodes'] # initial node coordinates
            self.ini_coord_faces = data['ini_faces'] # initial face coordinates
            self.nf = self.ini_coord_faces.shape[0]
            self.prev_disp = data['prev_nodes'] # previous total node displacement
            self.prev_face_disp = data['prev_faces'] # previous total face displacement
            self.dx = data['dx_nodes'] # latest time step node displacement
            self.face_dx = data['dx_faces'] # latest time step face displacement

        self.heat_flux = np.zeros((self.nf, 1)) # heat flux [W/m^2]

        # create input & output ModelParts
        self.model = Model()

        flag_nodes = False
        flag_faces = False
        for item in (self.settings['interface_input']):
            mp_name = item['model_part']
            if 'nodes' in mp_name:
                if flag_nodes:
                    raise ValueError('Only a single input model part for both nodes and faces is allowed for this Python solver.')
                else:
                    flag_nodes = True
                    self.input_nodes_mp_name = mp_name
                    self.model.create_model_part(self.input_nodes_mp_name, self.ini_coord_nodes[:, 0].flatten(),
                                                 self.ini_coord_nodes[:, 1].flatten(), np.zeros(self.nn),
                                                 np.arange(self.nn))
            elif 'faces' in mp_name:
                if flag_faces:
                    raise ValueError('Only a single input model part for both nodes and faces is allowed for this Python solver.')
                else:
                    flag_faces = True
                    self.input_faces_mp_name = mp_name
                    self.model.create_model_part(self.input_faces_mp_name, self.ini_coord_faces[:, 0].flatten(),
                                                 self.ini_coord_faces[:, 1].flatten(), np.zeros(self.nf), np.arange(self.nf))
            else:
                raise ValueError('Given model parts should contain "nodes" or "faces" indicating to which the variable is applied.')

        if not (flag_nodes and flag_faces):
            raise ValueError('One (and only one) "faces" model part (heat flux) and one "nodes" model part (displacement) is required as input for this Python solver.')

        if len(self.settings['interface_output']) == 1:
            self.output_mp_name = self.settings['interface_output'][0]['model_part']
            self.model.create_model_part(self.output_mp_name, self.ini_coord_nodes[:, 0].flatten(),
                                                   self.ini_coord_nodes[:, 1].flatten(), np.zeros(self.nn), np.arange(self.nn))
        else:
            raise ValueError('Only a single output model part for nodes (displacement) is allowed for this Python solver.')

        # input & output interfaces
        self.interface_input = Interface(self.settings['interface_input'], self.model)
        self.interface_input.set_variable_data(self.input_faces_mp_name, 'heat_flux', self.heat_flux)
        self.interface_input.set_variable_data(self.input_nodes_mp_name, 'displacement', self.prev_disp + self.dx)

        self.interface_output = Interface(self.settings['interface_output'], self.model)
        self.interface_output.set_variable_data(self.output_mp_name, 'displacement', self.prev_disp + self.dx)

        # internal models for face-to-node mapping
        self.internal_model = Model()
        self.internal_face_settings = [{"model_part": self.output_mp_name, "variables": ["1ts_disp", "prev_disp", "area"]}]
        self.internal_node_settings = [{"model_part": self.output_mp_name, "variables": ["displacement", "1ts_disp", "prev_disp"]}]

        self.internal_model.create_model_part(self.output_mp_name, self.ini_coord_faces[:, 0].flatten(),
                                     self.ini_coord_faces[:, 1].flatten(), np.zeros(self.nf), np.arange(self.nf))

        self.interface_internal_faces = Interface(self.internal_face_settings, self.internal_model)
        self.interface_internal_nodes = Interface(self.internal_node_settings, self.model)

        # create and initialize face to node (f2n) displacement mapper
        f2n_settings_dict = {"directions": ["x", "y"],
                             "check_bounding_box": False,
                             "projection_order": self.mapper_projection,
                             "domain": self.mapping_domain,
                             "limits": self.mapping_limits,
                             "mapping_conservative": self.conservative,
                             "mapper_relaxation": self.mapper_relaxation,
                             "mapper_iterations": self.mapper_iterations,
                             "volume_tolerance": self.vol_tolerance}

        f2n_settings = {"type": "mappers.interface", "settings": {"type": "mappers.linear_conservative",
                                                                  "settings": f2n_settings_dict}}
        self.mapper_f2n = tools.create_instance(f2n_settings)
        self.mapper_f2n.initialize(self.interface_internal_faces, self.interface_internal_nodes)

        # create and initialize node to face (n2f) displacement mapper
        n2f_settings = {"type": "mappers.interface", "settings": {"type": "mappers.linear",
                                                                  "settings": {"directions": ["x", "y"],
                                                                               "check_bounding_box": False}}}
        self.mapper_n2f = tools.create_instance(n2f_settings)
        self.mapper_n2f.initialize(self.interface_internal_nodes, self.interface_internal_faces)

    @tools.time_initialize
    def initialize(self):
        super().initialize()

    def initialize_solution_step(self):
        super().initialize_solution_step()
        self.timestep += 1

    @tools.time_solve_solution_step
    def solve_solution_step(self, interface_input):
        # process input interface data
        # store incoming variables
        self.interface_input.set_interface_data(interface_input.get_interface_data())

        self.heat_flux = self.interface_input.get_variable_data(self.input_faces_mp_name, 'heat_flux') # [W/m^2]
        self.heat_flux = -1 * self.heat_flux # negative heat flux for liquid domain is positive for solid domain

        # Shouldn't this be in initialize_solution_step??
        self.prev_disp = self.interface_input.get_variable_data(self.input_nodes_mp_name, 'displacement')
        self.interface_internal_nodes.set_variable_data(self.output_mp_name, 'prev_disp', self.prev_disp)

        print('\n')
        print(f'solve_solution_step - before n2f:')
        print(
            f'Norm interface_input.get_variable_data("displacement") = {np.linalg.norm(self.interface_input.get_variable_data(self.input_nodes_mp_name, "displacement"))}')
        print(
            f'Norm interface_internal_nodes.get_variable_data("prev_disp") = {np.linalg.norm(self.interface_internal_nodes.get_variable_data(self.output_mp_name, "prev_disp"))}')
        print(
            f'Norm interface_internal_nodes.get_variable_data("1ts_disp") = {np.linalg.norm(self.interface_internal_nodes.get_variable_data(self.output_mp_name, "1ts_disp"))}')
        print(
            f'Norm interface_internal_nodes.get_variable_data("displacement") = {np.linalg.norm(self.interface_internal_nodes.get_variable_data(self.output_mp_name, "displacement"))}')
        print(
            f'Norm interface_internal_faces.get_variable_data("prev_disp") = {np.linalg.norm(self.interface_internal_faces.get_variable_data(self.output_mp_name, "prev_disp"))}')
        print(
            f'Norm interface_internal_faces.get_variable_data("1ts_disp") = {np.linalg.norm(self.interface_internal_faces.get_variable_data(self.output_mp_name, "1ts_disp"))}')
        print(
            f'Norm interface_internal_faces.get_variable_data("area") = {np.linalg.norm(self.interface_internal_faces.get_variable_data(self.output_mp_name, "area"))}')

        # map previous node displacement to previous face displacement
        self.mapper_n2f.map_n2f(self.interface_internal_nodes, self.interface_internal_faces)

        print('\n')
        print(f'solve_solution_step - after n2f:')
        print(
            f'Norm interface_input.get_variable_data("displacement") = {np.linalg.norm(self.interface_input.get_variable_data(self.input_nodes_mp_name, "displacement"))}')
        print(
            f'Norm interface_internal_nodes.get_variable_data("prev_disp") = {np.linalg.norm(self.interface_internal_nodes.get_variable_data(self.output_mp_name, "prev_disp"))}')
        print(
            f'Norm interface_internal_nodes.get_variable_data("1ts_disp") = {np.linalg.norm(self.interface_internal_nodes.get_variable_data(self.output_mp_name, "1ts_disp"))}')
        print(
            f'Norm interface_internal_nodes.get_variable_data("displacement") = {np.linalg.norm(self.interface_internal_nodes.get_variable_data(self.output_mp_name, "displacement"))}')
        print(
            f'Norm interface_internal_faces.get_variable_data("prev_disp") = {np.linalg.norm(self.interface_internal_faces.get_variable_data(self.output_mp_name, "prev_disp"))}')
        print(
            f'Norm interface_internal_faces.get_variable_data("1ts_disp") = {np.linalg.norm(self.interface_internal_faces.get_variable_data(self.output_mp_name, "1ts_disp"))}')
        print(
            f'Norm interface_internal_faces.get_variable_data("area") = {np.linalg.norm(self.interface_internal_faces.get_variable_data(self.output_mp_name, "area"))}')

        disp_magn = (self.heat_flux * self.dt) / (self.rho * self.L) # Stefan condition
        self.area, normal_array = self.area_calc(self.ini_coord_nodes + self.prev_disp)

        self.face_dx = disp_magn * normal_array

        self.interface_internal_faces.set_variable_data(self.output_mp_name, '1ts_disp', self.face_dx)
        self.interface_internal_faces.set_variable_data(self.output_mp_name, 'area', self.area)

        print('\n')
        print(f'solve_solution_step - after calc, before f2n:')
        print(
            f'Norm interface_input.get_variable_data("displacement") = {np.linalg.norm(self.interface_input.get_variable_data(self.input_nodes_mp_name, "displacement"))}')
        print(
            f'Norm interface_internal_nodes.get_variable_data("prev_disp") = {np.linalg.norm(self.interface_internal_nodes.get_variable_data(self.output_mp_name, "prev_disp"))}')
        print(
            f'Norm interface_internal_nodes.get_variable_data("1ts_disp") = {np.linalg.norm(self.interface_internal_nodes.get_variable_data(self.output_mp_name, "1ts_disp"))}')
        print(
            f'Norm interface_internal_nodes.get_variable_data("displacement") = {np.linalg.norm(self.interface_internal_nodes.get_variable_data(self.output_mp_name, "displacement"))}')
        print(
            f'Norm interface_internal_faces.get_variable_data("prev_disp") = {np.linalg.norm(self.interface_internal_faces.get_variable_data(self.output_mp_name, "prev_disp"))}')
        print(
            f'Norm interface_internal_faces.get_variable_data("1ts_disp") = {np.linalg.norm(self.interface_internal_faces.get_variable_data(self.output_mp_name, "1ts_disp"))}')
        print(
            f'Norm interface_internal_faces.get_variable_data("area") = {np.linalg.norm(self.interface_internal_faces.get_variable_data(self.output_mp_name, "area"))}')

        # map face displacement to node displacement
        self.mapper_f2n.map_f2n(self.interface_internal_faces, self.interface_internal_nodes)
        
        self.dx = self.interface_internal_nodes.get_variable_data(self.output_mp_name, '1ts_disp')

        print('\n')
        print(f'solve_solution_step - after calc, after f2n:')
        print(
            f'Norm interface_input.get_variable_data("displacement") = {np.linalg.norm(self.interface_input.get_variable_data(self.input_nodes_mp_name, "displacement"))}')
        print(
            f'Norm interface_internal_nodes.get_variable_data("prev_disp") = {np.linalg.norm(self.interface_internal_nodes.get_variable_data(self.output_mp_name, "prev_disp"))}')
        print(
            f'Norm interface_internal_nodes.get_variable_data("1ts_disp") = {np.linalg.norm(self.interface_internal_nodes.get_variable_data(self.output_mp_name, "1ts_disp"))}')
        print(
            f'Norm interface_internal_nodes.get_variable_data("displacement") = {np.linalg.norm(self.interface_internal_nodes.get_variable_data(self.output_mp_name, "displacement"))}')
        print(
            f'Norm interface_internal_faces.get_variable_data("prev_disp") = {np.linalg.norm(self.interface_internal_faces.get_variable_data(self.output_mp_name, "prev_disp"))}')
        print(
            f'Norm interface_internal_faces.get_variable_data("1ts_disp") = {np.linalg.norm(self.interface_internal_faces.get_variable_data(self.output_mp_name, "1ts_disp"))}')
        print(
            f'Norm interface_internal_faces.get_variable_data("area") = {np.linalg.norm(self.interface_internal_faces.get_variable_data(self.output_mp_name, "area"))}')

        self.interface_output.set_variable_data(self.output_mp_name, 'displacement', self.prev_disp + self.dx)

        print('\n')
        print(f'solve_solution_step - output:')
        print(
            f'Norm interface_input.get_variable_data("displacement") = {np.linalg.norm(self.interface_input.get_variable_data(self.input_nodes_mp_name, "displacement"))}')
        print(
            f'Norm interface_output.get_variable_data("displacement") = {np.linalg.norm(self.interface_output.get_variable_data(self.output_mp_name, "displacement"))}')

        # output
        return self.interface_output

    def finalize_solution_step(self):
        super().finalize_solution_step()

    @tools.time_save
    def output_solution_step(self):
        super().output_solution_step()

        if self.timestep > 0 and self.save_restart != 0 and self.timestep % self.save_restart == 0:
            file_name = join(self.working_directory, f'case_timestep{self.timestep}.pickle')
            with open(file_name, 'wb') as file:
                pickle.dump({'ini_nodes': self.ini_coord_nodes,
                             'ini_faces': self.ini_coord_faces,
                             'prev_nodes': self.prev_disp,
                             'prev_faces': self.prev_face_disp,
                             'dx_nodes': self.dx,
                             'dx_faces': self.face_dx}, file)
            if self.save_restart < 0 and self.timestep + self.save_restart > self.timestep_start:
                try:
                    os.remove(join(self.working_directory, f'case_timestep{self.timestep + self.save_restart}.pickle'))
                except OSError:
                    pass

    def finalize(self):
        super().finalize()

    def area_calc(self, nodes):
        """
        Receives: array with ordered nodes (N, 3)
        Returns: array with ordered face areas (nf, 1) and their normals (nf, 3).
        """

        # 1. Define P1 (Start) and P2 (End) for all faces at once
        # Shift nodes by -1 so row 'i' aligns with row 'i+1'
        nodes_next = np.roll(nodes, -1, axis=0)

        if self.closed:
            # Closed: Use ALL points. Last point connects to First point (handled by roll).
            p1 = nodes
            p2 = nodes_next
        else:
            # Open: Exclude the last point from starts, and the wrapped point from ends.
            p1 = nodes[:-1]
            p2 = nodes_next[:-1]

        # 2. Vectorized Geometry Calculation
        # Calculate vector for every face (dx, dy, dz)
        diff = p2 - p1

        # Calculate lengths (Area) for all faces
        # We use norm of X and Y only (assuming 2D interface in 3D space)
        area = np.linalg.norm(diff[:, :2], axis=1)

        nf = len(area)
        normal_array = np.zeros((nf, 3))

        # 3. Calculate Normals (dy, -dx)
        # Create a mask to avoid division by zero
        valid_mask = area > 1e-15

        if np.any(valid_mask):
            dx = diff[valid_mask, 0]
            dy = diff[valid_mask, 1]
            lengths = area[valid_mask]

            # Standard normal: (dy, -dx) normalized
            normal_array[valid_mask, 0] = dy / lengths
            normal_array[valid_mask, 1] = -dx / lengths
            # Z stays 0

        # 4. Vectorized Direction Logic
        if isinstance(self.mov_dir, str):
            # Calculate Centroid and Face Centers
            centroid = np.mean(nodes, axis=0)
            face_centers = (p1 + p2) / 2
            vec_to_centroid = centroid - face_centers

            # Calculate Dot Product for all faces at once using einsum
            # (row-wise dot product of normal_array and vec_to_centroid)
            dot_products = np.einsum('ij,ij->i', normal_array, vec_to_centroid)

            if self.mov_dir == 'inward':
                # Flip if pointing away (dot < 0)
                flip_mask = dot_products < 0
                normal_array[flip_mask] *= -1

            elif self.mov_dir == 'outward':
                # Flip if pointing inward (dot > 0)
                flip_mask = dot_products > 0
                normal_array[flip_mask] *= -1

        else:
            # Vector-based direction (e.g., self.mov_dir = [0, 1])
            mov_dir_arr = np.array(self.mov_dir[:2])  # Ensure 2D

            # Matrix multiplication to get dot products for all faces
            dots = np.dot(normal_array[:, :2], mov_dir_arr)

            # Flip if opposed to movement direction
            flip_mask = dots < 0
            normal_array[flip_mask] *= -1

        return np.reshape(area, (nf, 1)), normal_array
