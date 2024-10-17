from coconut import data_structure
from coconut.coupling_components.solver_wrappers.solver_wrapper import SolverWrapper
from coconut import tools
from coconut.data_structure.variables import accepted_variables_pc

import os
from os.path import join
import glob
import subprocess
import multiprocessing
import numpy as np
import hashlib
from getpass import getuser
import shutil


def create(parameters):
    return SolverWrapperPCFluent(parameters)


class SolverWrapperPCFluent(SolverWrapper):
    # version specific parameters
    version = None  # Fluent product version, as from 2023R1 typically of the form 'xxxRx', set in subclass
    version_bis = None  # Fluent internal version, typically of the form 'x.x.0', set in subclass
    check_coupling_convergence_possible = False  # can solver check convergence after 1 iteration?

    @tools.time_initialize
    def __init__(self, parameters):
        super().__init__(parameters)

        if self.version is None or self.version_bis is None:
            raise NotImplementedError(
                'Base class method called, class variable version and version_bis need to be set in the derived class')

        # set parameters
        self.settings = parameters['settings']
        self.dir_cfd = join(os.getcwd(), self.settings['working_directory'])
        self.env = None  # environment in which correct version of Fluent software is available, set in sub-class
        self.coco_messages = tools.CocoMessages(self.dir_cfd)
        self.coco_messages.remove_all_messages()
        self.backup_fluent_log()
        self.dir_src = os.path.realpath(os.path.dirname(__file__))
        self.tmp_dir = os.environ.get('TMPDIR', '/tmp')  # dir for host-node communication
        self.tmp_dir_unique = os.path.join(self.tmp_dir, f'coconut_{getuser()}_{os.getpid()}_fluent')
        self.cores = self.settings['cores']
        self.hosts_file = self.settings.get('hosts_file')
        self.case_file = self.settings['case_file']
        self.data_file = self.case_file.replace('.cas', '.dat', 1)
        if not os.path.exists(os.path.join(self.dir_cfd, self.case_file)):
            raise FileNotFoundError(f'Case file {self.case_file} not found in working directory {self.dir_cfd}')
        elif not os.path.exists(os.path.join(self.dir_cfd, self.data_file)):
            raise FileNotFoundError(f'Data file {self.data_file} not found in working directory {self.dir_cfd}')
        self.mnpf = self.settings['max_nodes_per_face']
        self.dimensions = self.settings['dimensions']
        self.unsteady = self.settings['unsteady']
        self.multiphase = self.settings.get('multiphase', False)
        self.flow_iterations = self.settings['flow_iterations']
        self.delta_t = self.settings['delta_t']
        self.timestep_start = self.settings['timestep_start']
        self.timestep = self.timestep_start
        self.save_results = self.settings.get('save_results', 1)
        self.save_restart = self.settings['save_restart']
        self.iteration = None
        self.fluent_process = None
        self.thread_ids = {}  # thread IDs corresponding to thread names
        for thread_name in self.settings['thread_names']:
            self.thread_ids[thread_name] = None
        self.model_part_thread_ids = {}  # thread IDs corresponding to ModelParts
        self.dict_face_ids = {} # Dictionary of dictionaries containing a list of node ids corresponding to the hashed face ids
        self.model = None

        # Phase change specific settings
        self.pc_settings = self.settings['PC']
        self.ini_condition = self.pc_settings.get('ini_condition', None) # initial condition for outgoing variables
        self.latent_heat = self.pc_settings.get('latent', None)
        if self.latent_heat is None:
            raise ValueError('Latent heat should be provided in the json file in case of phase change problems')
        self.melt_temp = self.pc_settings.get('melt_temp', None)
        if self.melt_temp is None:
            raise ValueError('Melt temperature should be provided in the json file in case of phase change problems')
        self.moving_boundary = self.pc_settings.get('moving_boundary', True)
        self.mapping_domain = self.pc_settings.get('mapping_domain', None)
        self.mapping_limits = self.pc_settings.get('mapping_limits', None)

        self.output_ini_cond = {}
        self.output_variables = []
        f_n = 0
        f_f = 0
        for mp in self.settings['interface_output']:
            self.output_ini_cond[mp['model_part']] = 0
            if "nodes" in mp["model_part"]:
                if f_n == 0:
                    self.output_variables += mp['variables']  # NEW: list of output variables
                    f_n = 1
            if "faces" in mp["model_part"]:
                if f_f == 0:
                    self.output_variables += mp['variables']  # NEW: list of output variables
                    f_f = 1

        f_n = 0
        f_f = 0
        self.input_variables = []
        for mp in self.settings['interface_input']:
            if "nodes" in mp["model_part"]:
                if f_n == 0:
                    self.input_variables += mp['variables']  # NEW: list of input variables
                    f_n = 1
            if "faces" in mp["model_part"]:
                if f_f == 0:
                    self.input_variables += mp['variables']  # NEW: list of input variables
                    f_f = 1

        self.thermal_bc = None
        if "temperature" in self.input_variables:
            self.thermal_bc = "temperature"
        elif "heat_flux" in self.input_variables:
            self.thermal_bc = "heat_flux"

    @tools.time_initialize
    def initialize(self):
        super().initialize()

        # prepare Fluent journal
        journal = f'v{self.version}.jou'
        thread_names_str = ''
        for thread_name in self.thread_ids:
            thread_names_str += ' "' + thread_name + '"'
        if "pressure" in self.output_variables and "traction" in self.output_variables:
            force_balance = '#t'
        else:
            force_balance = '#f'
        thermal_bc = str(2)
        if self.thermal_bc == 'heat_flux':
            thermal_bc = str(1)
        elif self.thermal_bc == 'temperature':
            thermal_bc = str(0)
        with open(join(self.dir_src, journal)) as infile:
            with open(join(self.dir_cfd, journal), 'w') as outfile:
                for line in infile:
                    line = line.replace('|CASE|', join(self.dir_cfd, self.case_file))
                    line = line.replace('|THREAD_NAMES|', thread_names_str)
                    line = line.replace('|UNSTEADY|', '#t' if self.unsteady else '#f')
                    line = line.replace('|MULTIPHASE|', '#t' if self.multiphase else '#f')
                    line = line.replace('|MOVING_BOUNDARY|', '#t' if self.moving_boundary else '#f')
                    line = line.replace('|FORCE_BALANCE|', force_balance)
                    line = line.replace('|THERMAL_BC|', thermal_bc)
                    line = line.replace('|TEMP_STORED|', '#t' if "temperature" in self.output_variables else '#f')
                    line = line.replace('|FLOW_ITERATIONS|', str(self.flow_iterations))
                    line = line.replace('|DELTA_T|', str(self.delta_t))
                    line = line.replace('|TIMESTEP_START|', str(self.timestep_start))
                    line = line.replace('|END_OF_TIMESTEP_COMMANDS|', self.settings.get('end_of_timestep_commands','\n'))
                    line = line.replace('|END_OF_SETUP_COMMANDS|', self.pc_settings.get('end_of_setup_commands','\n'))
                    outfile.write(line)

        # prepare Fluent UDF
        udf = 'udf_thermal.c'
        with open(join(self.dir_src, udf)) as infile:
            with open(join(self.dir_cfd, udf), 'w') as outfile:
                for line in infile:
                    line = line.replace('|MAX_NODES_PER_FACE|', str(self.mnpf))
                    line = line.replace('|TMP_DIRECTORY_NAME|', self.tmp_dir_unique)
                    line = line.replace('|TIME_STEP_SIZE|', str(self.delta_t))
                    line = line.replace('|LATENT_HEAT|', str(self.latent_heat))
                    line = line.replace('|MELT_TEMP|', str(self.melt_temp))
                    line = line.replace('|UNSTEADY|', 'true' if self.unsteady else 'false')
                    line = line.replace('|FLUID|', 'false' if self.thermal_bc == "heat_flux" else 'true')
                    outfile.write(line)

        # check number of cores
        if self.hosts_file is not None:
            with open(join(self.dir_cfd, self.hosts_file)) as fp:
                max_cores = len(fp.readlines())
        else:
            max_cores = multiprocessing.cpu_count()
        if self.cores < 1 or self.cores > max_cores:
            warning = f'Number of cores incorrect, changed from {self.cores} to {max_cores}'
            if self.hosts_file is None:
                warning += f'\nAre you trying to run multinode?' \
                           f'\n\tAdd hosts file to working directory ({self.dir_cfd})' \
                           f'\n\tAdd "hosts_file" parameter with name of hosts file to Fluent parameters in JSON file' \
                           f'\nSee also https://pyfsi.github.io/coconut/fluent.html#running-multinode'
            else:
                warning += f'\nAre you trying to run multinode?' \
                           f'\n\tThe used hosts_file is {join(self.dir_cfd, self.hosts_file)}'
            tools.print_info(warning, layout='warning')
            self.cores = max_cores

        # check variables in parameter.json file
        for var in self.input_variables:
            if var not in accepted_variables_pc['in']:
                raise NameError("Only permitted input variables are displacement, temperature and heat flux")
        for var in self.output_variables:
            if var not in accepted_variables_pc['out']:
                raise NameError("Only permitted output variables are displacement, temperature, heat flux, pressure and traction")

        # start Fluent with journal
        log = join(self.dir_cfd, 'fluent.log')
        cmd1 = f'fluent -r{self.version_bis} {self.dimensions}ddp '
        cmd2 = f'-t{self.cores} -i {journal}'
        cmd3 = f' >> {log} 2>&1'

        if self.hosts_file is not None:
            cmd1 += f' -cnf={self.hosts_file} -ssh '
        if self.settings['fluent_gui']:
            cmd = cmd1 + cmd2 + cmd3
        else:
            cmd = cmd1 + '-gu ' + cmd2 + cmd3
        self.fluent_process = subprocess.Popen(cmd, executable='/bin/bash',
                                               shell=True, cwd=self.dir_cfd, env=self.env)

        # get general simulation info from  fluent.log and report.sum
        self.coco_messages.wait_message('case_info_exported')

        with open(log, 'r') as file:
            for line in file:
                if 'File has wrong dimension' in line:
                    raise ValueError('Dimension in JSON does not match Fluent case')

        report = join(self.dir_cfd, 'report.sum')
        check = 0
        with open(report, 'r') as file:
            for line in file:
                if 'Model' in line and 'Settings' in line:
                    check = 1
                elif check == 1 and 'Space' in line:
                    if str(self.dimensions) not in line:
                        if not (self.dimensions == 2 and 'Axisymmetric' in line):
                            raise ValueError(f'Dimension in JSON does not match Fluent')
                    check = 2
                elif check == 2 and 'Time' in line:
                    if 'Steady' in line and self.unsteady:
                        raise ValueError('Unsteady in JSON does not match steady Fluent')
                    elif 'Unsteady' in line and not self.unsteady:
                        raise ValueError('Steady in JSON does not match unsteady Fluent')
                    check = 3
                elif check == 3 and 'Equation' in line and 'Solved' in line:
                    check = 4
                elif check == 4:
                    if 'Volume Fraction' in line and 'yes' in line:
                        if not self.multiphase:
                            raise ValueError('Singlephase in JSON does not match multiphase Fluent')
                        break
                    elif 'Numerics' in line:
                        if self.multiphase:
                            raise ValueError('Multiphase in JSON does not match singlephase Fluent')
                        break

        if os.path.isfile(join(self.dir_cfd, 'log')):
            os.unlink(join(self.dir_cfd, 'log'))  # delete log file (fluent.log is sufficient)

        # get surface thread ID's from report.sum and write them to bcs.txt
        check = 0
        names_found = []
        with open(report, 'r') as file:
            for line in file:
                if check == 3 and line.islower():
                    line_list = line.strip().split()
                    if len(line_list) == 3:
                        name, thread_id, _ = line_list
                    elif len(line_list) == 4:
                        name, _, thread_id, _ = line_list
                    else:
                        raise ValueError(f'Format of {report} not recognized')
                    if name in self.thread_ids and name not in names_found:
                        self.thread_ids[name] = thread_id
                        names_found.append(name)
                if check == 3 and not line.islower():
                    break
                if check == 2:  # skip 1 line
                    check = 3
                if 'name' in line and check == 1:
                    check = 2
                if 'Boundary Conditions' in line:
                    check = 1
        with open(join(self.dir_cfd, 'bcs.txt'), 'w') as file:
            file.write(f'{len(names_found)}\n')
            for name, id in self.thread_ids.items():
                file.write(f'{name} {id}\n')
        self.coco_messages.send_message('thread_ids_written_to_file')

        # remove "report.sum" because the batch options to overwrite report files and case files conflict in some
        # versions of Fluent (2023R1)
        os.unlink(report)

        # import node and face information if no restart
        if self.timestep_start == 0:
            self.coco_messages.wait_message('nodes_and_faces_stored')

        # create Model used for coupling
        self.model = data_structure.Model()

        # create Model used internally for face displacements
        if "displacement" in self.output_variables:
            self.internal_model = data_structure.Model()
            self.internal_settings = []

        # create input ModelParts (nodes for displacement, faces for temperature or heat flux)
        for item in (self.settings['interface_input']):
            mp_name = item['model_part']

            # get face thread ID that corresponds to ModelPart
            for thread_name in self.thread_ids:
                if thread_name in mp_name:
                    self.model_part_thread_ids[mp_name] = self.thread_ids[thread_name]
            if mp_name not in self.model_part_thread_ids:
                raise AttributeError('Could not find thread name corresponding ' +
                                     f'to ModelPart {mp_name}')

            # read in datafile
            thread_id = self.model_part_thread_ids[mp_name]
            if "nodes" in mp_name:
                file_name = join(self.dir_cfd, f'nodes_timestep0_thread{thread_id}.dat')
                data = np.loadtxt(file_name, skiprows=1)
                if data.shape[1] != self.dimensions + 1:
                    raise ValueError('Given dimension does not match coordinates')
            elif "faces" in mp_name:
                file_name = join(self.dir_cfd, f'faces_timestep0_thread{thread_id}.dat')
                data = np.loadtxt(file_name, skiprows=1, ndmin=2)
                if data.shape[1] != self.dimensions + self.mnpf:
                    raise ValueError(f'Given dimension does not match coordinates')

            # get node or face coordinates and ids
            coords_tmp = np.zeros((data.shape[0], 3)) * 0.
            if "nodes" in mp_name:
                coords_tmp[:, :self.dimensions] = data[:, :-1]  # add column z if 2D
                ids_tmp = data[:, -1].astype(int)  # array is flattened
            elif "faces" in mp_name:
                coords_tmp[:, :self.dimensions] = data[:, :-self.mnpf]  # add column z if 2D
                ids_tmp, self.dict_face_ids[mp_name] = self.get_unique_face_ids(data[:, -self.mnpf:])

            # sort and remove doubles
            args = np.unique(ids_tmp, return_index=True)[1].tolist()
            x0 = coords_tmp[args, 0]
            y0 = coords_tmp[args, 1]
            z0 = coords_tmp[args, 2]
            ids = ids_tmp[args]

            # create ModelPart
            self.model.create_model_part(mp_name, x0, y0, z0, ids)

        # create output ModelParts (nodes or faces)
        for item in self.settings['interface_output']:
            mp_name = item['model_part']

            # get face thread ID that corresponds to ModelPart
            for thread_name in self.thread_ids:
                if thread_name in mp_name:
                    self.model_part_thread_ids[mp_name] = self.thread_ids[thread_name]
            if mp_name not in self.model_part_thread_ids:
                raise AttributeError('Could not find thread name corresponding ' +
                                     f'to ModelPart {mp_name}')

            # read in datafile
            thread_id = self.model_part_thread_ids[mp_name]
            # face coordinates are necessary for both face and nodes variables,
            # because displacement in nodes must first be processed at face level
            file_name = join(self.dir_cfd, f'faces_timestep0_thread{thread_id}.dat')
            data = np.loadtxt(file_name, skiprows=1, ndmin=2)
            if data.shape[1] != self.dimensions + self.mnpf:
                raise ValueError(f'Given dimension does not match coordinates')
            if "nodes" in mp_name:
                file_name = join(self.dir_cfd, f'nodes_timestep0_thread{thread_id}.dat')
                data_nodes = np.loadtxt(file_name, skiprows=1, ndmin=2)
                if data_nodes.shape[1] != self.dimensions + 1:
                    raise ValueError('Given dimension does not match coordinates')

            # get node or face coordinates and ids
            coords_tmp = np.zeros((data.shape[0], 3)) * 0.
            coords_tmp[:, :self.dimensions] = data[:, :-self.mnpf]  # add column z if 2D
            ids_tmp, self.dict_face_ids[mp_name] = self.get_unique_face_ids(data[:, -self.mnpf:])
            if "nodes" in mp_name:
                coords_nodes_tmp = np.zeros((data_nodes.shape[0], 3)) * 0.
                coords_nodes_tmp[:, :self.dimensions] = data_nodes[:, :-1]  # add column z if 2D
                ids_nodes_tmp = data_nodes[:, -1].astype(int)  # array is flattened

            # sort and remove doubles
            args = np.unique(ids_tmp, return_index=True)[1].tolist()
            x0 = coords_tmp[args, 0]
            y0 = coords_tmp[args, 1]
            z0 = coords_tmp[args, 2]
            ids = ids_tmp[args]
            if "nodes" in mp_name:
                args_nodes = np.unique(ids_nodes_tmp, return_index=True)[1].tolist()
                x0_n = coords_nodes_tmp[args_nodes, 0]
                y0_n = coords_nodes_tmp[args_nodes, 1]
                z0_n = coords_nodes_tmp[args_nodes, 2]
                ids_n = ids_nodes_tmp[args_nodes]

            # create ModelPart
            if "faces" in mp_name:
                self.model.create_model_part(mp_name, x0, y0, z0, ids)
            elif "nodes" in mp_name:
                self.model.create_model_part(mp_name, x0_n, y0_n, z0_n, ids_n)
                # mp in internal_model has "nodes" in name, but has coordinates and ids of faces
                self.internal_settings.append({"model_part": mp_name, "variables": ["displacement"]})
                self.internal_model.create_model_part(mp_name, x0, y0, z0, ids)

            # create initial conditions at output interface
            if self.ini_condition is not None:
                if "faces" in mp_name:
                    self.output_ini_cond[mp_name] = np.ones((data.shape[0], 1))*self.ini_condition

        # create Interfaces
        self.interface_input = data_structure.Interface(self.settings['interface_input'], self.model)
        self.interface_output = data_structure.Interface(self.settings['interface_output'], self.model)
        if "displacement" in self.output_variables:
            self.interface_internal = data_structure.Interface(self.internal_settings, self.internal_model)

        # set initial conditions at output interface
        for key in self.output_ini_cond:
            if "faces" in key:
                for var in self.output_variables:
                    if var == 'heat_flux' or var == 'temperature':
                        if self.ini_condition is None:
                            raise ValueError('Initial condition must be defined in JSON for temperature or heat flux')
                        else:
                            self.interface_output.set_variable_data(key, var, self.output_ini_cond[key])

        if "displacement" in self.output_variables:
            # create, initialize and query face to node (f2n) displacement mapper
            mapper_settings = {"type": "mappers.interface", "settings": {"type": "mappers.linear",
                                                                         "settings": {"directions": ["x", "y"],
                                                                                      "check_bounding_box": False,
                                                                                      "domain": self.mapping_domain,
                                                                                      "limits": self.mapping_limits}}}
            self.mapper_f2n = tools.create_instance(mapper_settings)
            self.mapper_f2n.initialize(self.interface_internal, self.interface_output)
            # create, initialize and query node to face (n2f) displacement mapper
            self.mapper_n2f = tools.create_instance(mapper_settings)
            self.mapper_n2f.initialize(self.interface_output, self.interface_internal)

    def initialize_solution_step(self):
        super().initialize_solution_step()

        self.iteration = 0
        self.timestep += 1

        if "displacement" in self.output_variables:
            # map output node displacement to face displacement
            self.mapper_n2f.map_n2f(self.interface_output, self.interface_internal)

            # write interface output data at new time step
            self.write_output_to_file("displacement")

        self.coco_messages.send_message('next')
        self.coco_messages.wait_message('next_ready')

    @tools.time_solve_solution_step
    def solve_solution_step(self, interface_input):
        self.iteration += 1

        # process input interface data
        # store incoming variables
        self.interface_input.set_interface_data(interface_input.get_interface_data())

        for var in self.input_variables:
            # write interface input data
            self.write_input_to_file(var)

            # copy input data for debugging
            if self.debug:
                for dct in self.interface_input.parameters:
                    mp_name = dct['model_part']
                    thread_id = self.model_part_thread_ids[mp_name]
                    src = accepted_variables_pc['in'][var][0] + f'_timestep{self.timestep}_thread{thread_id}.dat'
                    dst = accepted_variables_pc['in'][var][0] + f'_timestep{self.timestep}_thread{thread_id}_Iter{self.iteration}.dat'
                    cmd = f'cp {join(self.dir_cfd, src)} {join(self.dir_cfd, dst)}'
                    os.system(cmd)

        # let Fluent run, wait for data
        self.coco_messages.send_message('continue')
        self.coco_messages.wait_message('continue_ready')

        # process output interface data
        for dct in self.interface_output.parameters:
            mp_name = dct['model_part']
            thread_id = self.model_part_thread_ids[mp_name]

            # read in datafile
            for var in dct['variables']:
                prefix = accepted_variables_pc['out'][var][0]
                # Avoid repeat of commands in case variables are stored in the same file
                if prefix != 'repeat':
                    data = self.read_output_file(prefix, mp_name, thread_id)
                    if prefix == "displacement" and self.thermal_bc != "heat_flux":
                        # Displacement is directly known at nodes in liquid solver
                        req_dim = accepted_variables_pc['out'][var][1] + 1
                    elif accepted_variables_pc['out'][var][1] == 0:
                        req_dim = self.dimensions + 1 + self.mnpf
                    else:
                        req_dim = min([self.dimensions, accepted_variables_pc['out'][var][1]]) + self.mnpf
                    if data.shape[1] != req_dim:
                        raise ValueError('Given dimension does not match coordinates')

                    if accepted_variables_pc['out'][var][1] == 0:
                        # get face coordinates and ids
                        vector_tmp = np.zeros((data.shape[0], 3)) * 0.
                        vector_tmp[:, :self.dimensions] = data[:, :-1 - self.mnpf]
                        scalar_tmp = data[:, self.dimensions].reshape(-1, 1)
                        ids_tmp, _ = self.get_unique_face_ids(data[:, -self.mnpf:])

                        # sort and remove doubles
                        args = np.unique(ids_tmp, return_index=True)[1].tolist()
                        vector = vector_tmp[args, :]
                        scalar = scalar_tmp[args]
                        ids = ids_tmp[args]

                        # store vector and scalar values in interface
                        model_part = self.model.get_model_part(mp_name)
                        if ids.size != model_part.size:
                            raise ValueError('Size of data does not match size of ModelPart')
                        if not np.all(ids == model_part.id):
                            raise ValueError('IDs of data do not match ModelPart IDs')

                        if var == 'pressure':
                            self.interface_output.set_variable_data(mp_name, 'traction', vector)
                            self.interface_output.set_variable_data(mp_name, var, scalar)
                        elif var == 'traction':
                            self.interface_output.set_variable_data(mp_name, var, vector)
                            self.interface_output.set_variable_data(mp_name, 'pressure', scalar)
                        elif var == 'heat_flux':
                            self.interface_output.set_variable_data(mp_name, var, scalar)

                    if accepted_variables_pc['out'][var][1] == 3:
                        # get face coordinates and ids
                        vector_tmp = np.zeros((data.shape[0], 3)) * 0.
                        vector_tmp[:, :self.dimensions] = data[:, : -self.mnpf]
                        ids_tmp, _ = self.get_unique_face_ids(data[:, -self.mnpf:])

                        # sort and remove doubles
                        args = np.unique(ids_tmp, return_index=True)[1].tolist()
                        vector = vector_tmp[args, :]
                        ids = ids_tmp[args]

                        if var == 'displacement' and self.thermal_bc == "heat_flux":
                            if "nodes" not in mp_name:
                                raise ValueError('Model part must be node-based for the displacement variable')

                            displacement = self.interface_internal.get_variable_data(mp_name, 'displacement')
                            displacement_new = displacement + vector

                            model_part = self.internal_model.get_model_part(mp_name)
                            if ids.size != model_part.size:
                                raise ValueError('Size of data does not match size of ModelPart')
                            if not np.all(ids == model_part.id):
                                raise ValueError('IDs of data do not match ModelPart IDs')

                            self.interface_internal.set_variable_data(mp_name, var, displacement_new)

                            # map face displacement to node displacement
                            self.mapper_f2n.map_f2n(self.interface_internal, self.interface_output)

                            # write nodes_update file for move nodes udf in solid domain
                            self.write_output_to_file(var)

                            # reset face displacement to converged position of previous time step (avoid accumulation)
                            self.interface_internal.set_variable_data(mp_name, var, displacement)

                        else:
                            # store vector values in interface
                            model_part = self.model.get_model_part(mp_name)
                            if ids.size != model_part.size:
                                raise ValueError('Size of data does not match size of ModelPart')
                            if not np.all(ids == model_part.id):
                                raise ValueError('IDs of data do not match ModelPart IDs')

                            self.interface_output.set_variable_data(mp_name, var, vector)

                    if accepted_variables_pc['out'][var][1] == 1:
                        # get face coordinates and ids
                        scalar_tmp = data[:, 0].reshape(-1, 1)
                        ids_tmp, _ = self.get_unique_face_ids(data[:, -self.mnpf:])

                        # sort and remove doubles
                        args = np.unique(ids_tmp, return_index=True)[1].tolist()
                        scalar = scalar_tmp[args]
                        ids = ids_tmp[args]

                        # store scalar values in interface
                        model_part = self.model.get_model_part(mp_name)
                        if ids.size != model_part.size:
                            raise ValueError('Size of data does not match size of ModelPart')
                        if not np.all(ids == model_part.id):
                            raise ValueError('IDs of data do not match ModelPart IDs')

                        self.interface_output.set_variable_data(mp_name, var, scalar)

        # return interface_output object
        return self.interface_output

    def finalize_solution_step(self):
        super().finalize_solution_step()

        # save if required
        if (self.save_results != 0 and self.timestep % self.save_results == 0) \
                or (self.save_restart != 0 and self.timestep % self.save_restart == 0):
            self.coco_messages.send_message('save')
            self.coco_messages.wait_message('save_ready')

        # remove unnecessary files
        if self.timestep - 1 > self.timestep_start:
            self.remove_dat_files(self.timestep - 1)
            if self.save_restart < 0 and self.timestep + self.save_restart > self.timestep_start and \
                    self.timestep % self.save_restart == 0 \
                    and (self.save_results == 0 or (self.timestep + self.save_restart) % self.save_results != 0):
                # new restart file is written (self.timestep % self.save_restart ==0),
                # so previous one (at self.timestep + self.save_restart) can be deleted if:
                # - save_restart is negative
                # - files from a previous calculation are not touched
                # - files are not kept for save_results

                for extension in ('cas.h5', 'cas', 'dat.h5', 'dat'):
                    try:
                        os.remove(join(self.dir_cfd, f'case_timestep{self.timestep + self.save_restart}.{extension}'))
                    except OSError:
                        continue

    def finalize(self):
        super().finalize()
        shutil.rmtree(self.tmp_dir_unique, ignore_errors=True)
        self.coco_messages.send_message('stop')
        self.coco_messages.wait_message('stop_ready')
        self.fluent_process.wait()

        # remove unnecessary files
        self.remove_dat_files(self.timestep)

        # delete .trn files
        if not self.debug:
            for path in glob.glob(join(self.dir_cfd, '*.trn')):
                os.remove(path)

    def remove_dat_files(self, timestep):
        if not self.debug:
            for thread_id in self.thread_ids.values():
                for var in self.output_variables:
                    if accepted_variables_pc['out'][var][0] != 'repeat':
                        try:
                            os.remove(join(self.dir_cfd, accepted_variables_pc['out'][var][0] + f'_timestep{timestep}_thread{thread_id}.dat'))
                        except OSError:
                            pass
                    if var == 'displacement':
                        # solid zone outputs displacement and creates input file for itself at the same time
                        try:
                            os.remove(join(self.dir_cfd, accepted_variables_pc['in'][var][0] + f'_timestep{timestep}_thread{thread_id}.dat'))
                        except OSError:
                            pass
                for var in self.input_variables:
                    try:
                        os.remove(join(self.dir_cfd, accepted_variables_pc['in'][var][0] + f'_timestep{timestep}_thread{thread_id}.dat'))
                    except OSError:
                        pass

    def check_software(self):
        # Fluent version: see set_fluent_version
        result = subprocess.run(['fluent', '-r'], stdout=subprocess.PIPE, env=self.env)
        if self.version_bis not in str(result.stdout):
            raise RuntimeError(f'ANSYS Fluent version {self.version} ({self.version_bis}) is required. Check if '
                               f'the solver load commands for the "machine_name" are correct in solver_modules.py.')

    # noinspection PyMethodMayBeStatic
    def get_unique_face_ids(self, data):
        """
        Construct unique face IDs based on the face's node IDs.

        Parameter data contains a 2D ndarray of node IDs.
        Each row corresponds to the unique node IDs corresponding
        to a certain face, supplemented with -1-values.
        The row is sorted, the -1-values are removed, and then a
        string is made by adding the unique node IDs together.
            e.g. for a row [5, 9, 7, -1, -1]
                 the face ID is "5-7-9"

        # *** NEW: as we need integer IDs, the string is hashed and shortened
        # *** NEW: store the unique string alongside the hash value to enable reconstruction
        """
        data = data.astype(int)
        # ids = np.zeros(data.shape[0], dtype='U256')  # array is flattened
        ids = np.zeros(data.shape[0], dtype=int)
        face_ids = {}  # Dictionary to store face IDs and corresponding unique strings
        for i in range(ids.size):
            tmp = np.unique(data[i, :])
            if tmp[0] == -1:
                tmp = tmp[1:]
            unique_string = '-'.join(tuple(tmp.astype(str)))
            hash_id = hashlib.sha1(str.encode(unique_string))
            ids[i] = int(hash_id.hexdigest(), 16) % (10 ** 16)
            face_ids[ids[i]] = unique_string
        return ids, face_ids

    def reverse_face_ids(self, id, mp_name):
        str_ids = self.dict_face_ids[mp_name][id]
        int_strings = str_ids.split("-")
        node_ids = [int(x) for x in int_strings]
        return node_ids

    def write_input_to_file(self, var):
        if var == 'displacement':
            for dct in self.interface_input.parameters:
                mp_name = dct['model_part']
                if 'nodes' in mp_name:
                    thread_id = self.model_part_thread_ids[mp_name]
                    model_part = self.model.get_model_part(mp_name)
                    displacement = self.interface_input.get_variable_data(mp_name, 'displacement')
                    x = model_part.x0 + displacement[:, 0]
                    y = model_part.y0 + displacement[:, 1]
                    z = model_part.z0 + displacement[:, 2]
                    if self.dimensions == 2:
                        data = np.rec.fromarrays([x, y, model_part.id])
                        fmt = '%27.17e%27.17e%27d'
                    else:
                        data = np.rec.fromarrays([x, y, z, model_part.id])
                        fmt = '%27.17e%27.17e%27.17e%27d'
                    tmp = f'nodes_update_timestep{self.timestep}_thread{thread_id}.dat'
                    file_name = join(self.dir_cfd, tmp)
                    np.savetxt(file_name, data, fmt=fmt, header=f'{model_part.size}', comments='')
        if var == 'temperature':
            for dct in self.interface_input.parameters:
                mp_name = dct['model_part']
                if 'faces' in mp_name:
                    thread_id = self.model_part_thread_ids[mp_name]
                    model_part = self.model.get_model_part(mp_name)
                    T = self.interface_input.get_variable_data(mp_name, 'temperature')
                    face_nodeIDs = np.zeros((np.size(model_part.id), self.mnpf))
                    for index, id in enumerate(model_part.id):
                        list = self.reverse_face_ids(id, mp_name)
                        if len(list) == self.mnpf:
                            face_nodeIDs[index,:] = np.array(list)
                        else:
                            raise ValueError("Nr. of nodes returned from reverse_face_ids function does not correspond to mnpf")
                    prof = np.append(T, face_nodeIDs, axis=1)
                    fmt = '%27.17e'
                    for i in range(self.mnpf):
                        fmt += '%27d'
                    tmp = f'temperature_timestep{self.timestep}_thread{thread_id}.dat'
                    file_name = join(self.dir_cfd, tmp)
                    np.savetxt(file_name, prof, fmt=fmt, header='temperature unique-ids', comments='')
        if var == 'heat_flux':
            for dct in self.interface_input.parameters:
                mp_name = dct['model_part']
                if 'faces' in mp_name:
                    thread_id = self.model_part_thread_ids[mp_name]
                    model_part = self.model.get_model_part(mp_name)
                    q_tmp = self.interface_input.get_variable_data(mp_name, 'heat_flux')
                    q = np.append(np.zeros((np.size(q_tmp), self.dimensions)), q_tmp, axis=1)
                    face_nodeIDs = np.zeros((np.size(model_part.id), self.mnpf))
                    for index, id in enumerate(model_part.id):
                        list = self.reverse_face_ids(id, mp_name)
                        if len(list) == self.mnpf:
                            face_nodeIDs[index,:] = np.array(list)
                        else:
                            raise ValueError("Nr. of nodes returned from reverse_face_ids function does not correspond to mnpf")
                    prof = np.append(q, face_nodeIDs, axis=1)
                    if self.dimensions == 2:
                        fmt = '%27.17e%27.17e%27.17e'
                        header = 'x-flux y-flux flux-normal unique-ids'
                    else:
                        fmt = '%27.17e%27.17e%27.17e%27.17e'
                        header = 'x-flux y-flux z-flux flux-normal unique-ids'
                    for i in range(self.mnpf):
                        fmt += '%27d'
                    tmp = f'heat_flux_timestep{self.timestep}_thread{thread_id}.dat'
                    file_name = join(self.dir_cfd, tmp)
                    np.savetxt(file_name, prof, fmt=fmt, header=header, comments='')

    def write_output_to_file(self, var):
        if var == 'displacement':
            for dct in self.interface_output.parameters:
                mp_name = dct['model_part']
                if 'nodes' in mp_name:
                    thread_id = self.model_part_thread_ids[mp_name]
                    model_part = self.model.get_model_part(mp_name)
                    displacement = self.interface_output.get_variable_data(mp_name, 'displacement')
                    x = model_part.x0 + displacement[:, 0]
                    y = model_part.y0 + displacement[:, 1]
                    z = model_part.z0 + displacement[:, 2]
                    if self.dimensions == 2:
                        data = np.rec.fromarrays([x, y, model_part.id])
                        fmt = '%27.17e%27.17e%27d'
                    else:
                        data = np.rec.fromarrays([x, y, z, model_part.id])
                        fmt = '%27.17e%27.17e%27.17e%27d'
                    tmp = f'nodes_update_timestep{self.timestep}_thread{thread_id}.dat'
                    file_name = join(self.dir_cfd, tmp)
                    np.savetxt(file_name, data, fmt=fmt, header=f'{model_part.size}', comments='')

    def read_output_file(self, prefix, mp_name, thread_id=None):
        tmp = prefix + f'_timestep{self.timestep}_thread{thread_id}.dat'
        file_name = join(self.dir_cfd, tmp)
        data = np.loadtxt(file_name, skiprows=1, ndmin=2)
        # copy output data for debugging
        if self.debug:
            dst = prefix + f'_timestep{self.timestep}_thread{thread_id}_it{self.iteration}.dat'
            cmd = f'cp {file_name} {join(self.dir_cfd, dst)}'
            os.system(cmd)
        return data

    def get_coordinates(self):
        """  # TODO: rewrite this + include input ModelParts for faces (only used in Fluent solver wrapper tests atm)
        This function can be used e.g. for debugging or testing.
        It returns a dict that contains keys for the ModelParts
        on the two Interfaces.
        These keys give other dicts that have keys 'ids' and
        'coords'. These refer to ndarrays with respectively the
        ids of all Nodes and the current coordinates of those
        Nodes in Fluent (i.e. of the deformed geometry).
        """

        # make Fluent store coordinates and ids
        self.coco_messages.send_message('store_grid')
        self.coco_messages.wait_message('store_grid_ready')

        coord_data = {}

        # get ids and coordinates for input ModelParts (nodes)
        for dct in self.interface_input.parameters:
            mp_name = dct['model_part']
            coord_data[mp_name] = {}
            thread_id = self.model_part_thread_ids[mp_name]

            # read in datafile
            tmp = f'nodes_timestep{self.timestep}_thread{thread_id}.dat'
            data = np.loadtxt(join(self.dir_cfd, tmp), skiprows=1)
            if data.shape[1] != self.dimensions + 1:
                raise ValueError('Given dimension does not match coordinates')

            # get node coordinates and ids
            coords_tmp = np.zeros((data.shape[0], 3)) * 0.
            coords_tmp[:, :self.dimensions] = data[:, :-1]  # add column z if 2D
            ids_tmp = data[:, -1].astype(int)  # array is flattened

            # sort and remove doubles
            args = np.unique(ids_tmp, return_index=True)[1].tolist()
            coord_data[mp_name]['ids'] = ids_tmp[args]
            coord_data[mp_name]['coords'] = coords_tmp[args, :]

        # update coordinates for output ModelParts (faces)
        for dct in self.interface_output.parameters:
            mp_name = dct['model_part']
            coord_data[mp_name] = {}
            thread_id = self.model_part_thread_ids[mp_name]

            # read in datafile
            tmp = f'faces_timestep{self.timestep}_thread{thread_id}.dat'
            data = np.loadtxt(join(self.dir_cfd, tmp), skiprows=1, ndmin=2)
            if data.shape[1] != self.dimensions + self.mnpf:
                raise ValueError(f'Given dimension does not match coordinates')

            # get face coordinates and ids
            coords_tmp = np.zeros((data.shape[0], 3)) * 0.
            coords_tmp[:, :self.dimensions] = data[:, :-self.mnpf]  # add column z if 2D
            ids_tmp, _ = self.get_unique_face_ids(data[:, -self.mnpf:])

            # sort and remove doubles
            args = np.unique(ids_tmp, return_index=True)[1].tolist()
            coord_data[mp_name]['ids'] = ids_tmp[args]
            coord_data[mp_name]['coords'] = coords_tmp[args, :]

        return coord_data

    def backup_fluent_log(self):
        file = join(self.dir_cfd, 'fluent.log')
        file_backup = join(self.dir_cfd, 'fluent_backup.log')
        if os.path.isfile(file_backup):
            os.remove(file_backup)
        if os.path.isfile(file):
            os.rename(file, file_backup)
