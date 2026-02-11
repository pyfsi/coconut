from coconut import data_structure
from coconut.coupling_components.solver_wrappers.solver_wrapper import SolverWrapper
from coconut import tools
from coconut.data_structure.variables import accepted_variables_pc_liquid_rb

import os
from os.path import join
import glob
import subprocess
import multiprocessing
import numpy as np
import math as m
import pandas as pd
from scipy.spatial import cKDTree
import hashlib
from getpass import getuser
import shutil
import pickle
import csv

def create(parameters):
    return SolverWrapperPCFluentLiquidRB(parameters)


class SolverWrapperPCFluentLiquidRB(SolverWrapper):
    # version specific parameters
    version = None  # Fluent product version, as from 2023R1 typically of the form 'xxxRx', set in subclass
    version_bis = None  # Fluent internal version, typically of the form 'x.x.0', set in subclass
    check_coupling_convergence_possible = False  # can solver check convergence after 1 iteration?

    # define input and output variables
    accepted_in_var = list(accepted_variables_pc_liquid_rb['in'].keys())
    accepted_out_var = list(accepted_variables_pc_liquid_rb['out'].keys())

    @tools.time_initialize
    def __init__(self, parameters):
        super().__init__(parameters)

        if self.version is None or self.version_bis is None:
            raise NotImplementedError(
                'Base class method called, class variable version and version_bis need to be set in the derived class')

        # set parameters
        self.settings = parameters['settings']
        self.dir_cfd = join(os.getcwd(), self.settings['working_directory'])
        self.env = None  # environment in which correct version of Fluent software is available, set in subclass
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
        self.restart = self.timestep_start != 0  # true if restart
        self.save_results = self.settings.get('save_results', 1)
        self.save_restart = self.settings['save_restart']
        self.iteration = None
        self.rb_iter = None
        self.fluent_process = None
        self.thread_ids = {}  # thread IDs corresponding to thread names
        for thread_name in self.settings['thread_names']:
            self.thread_ids[thread_name] = None
        self.model_part_thread_ids = {}  # thread IDs corresponding to ModelParts
        self.dict_face_ids = {} # Dictionary of dictionaries containing a list of node ids corresponding to the hashed face ids
        self.model = None
        self.outer_inner_surf_dict = self.settings['overset_boundaries']
        self.overset_thread_ids = {}
        self.define_os = True
        if len(self.outer_inner_surf_dict) == 0:
            self.define_os = False
        else:
            for thread_name in self.outer_inner_surf_dict.values():
                self.overset_thread_ids[thread_name] = None

        # Input variables
        f_n = 0
        f_f = 0
        self.input_variables = []
        check_in_out_disp = []
        for mp in self.settings['interface_input']:
            if "nodes" in mp["model_part"]:
                if f_n == 0:
                    self.input_variables += mp['variables']
                    f_n = 1
                if mp['variables'] == 'displacement':
                    check_in_out_disp.append(mp["model_part"])
            if "faces" in mp["model_part"]:
                if f_f == 0:
                    self.input_variables += mp['variables']
                    f_f = 1

        # Output variables
        self.output_ini_cond = {}
        self.output_variables = []
        f_n = 0
        f_f = 0
        for mp in self.settings['interface_output']:
            self.output_ini_cond[mp['model_part']] = 0
            if "nodes" in mp["model_part"]:
                if f_n == 0:
                    self.output_variables += mp['variables']
                    f_n = 1
                if mp['variables'] == 'displacement':
                    mp_in = mp["model_part"].replace('out', 'in')
                    if mp_in in check_in_out_disp:
                        check_in_out_disp.remove(mp_in)
                    else:
                        raise ValueError(f'Model part with output displacement {mp["model_part"]} lacks corresponding input model part: {mp_in}.')
            if "faces" in mp["model_part"]:
                if f_f == 0:
                    self.output_variables += mp['variables']
                    f_f = 1

        # Check if input model parts with displacement variable have no corresponding output model parts.
        # This is necessary for the interface update in the solid solver.
        if check_in_out_disp:
            absent_out_mp = [str(mp).replace('in', 'out') for mp in check_in_out_disp]
            raise ValueError(f'Model parts with input displacement {check_in_out_disp} lack corresponding output model parts to return the displacement: {absent_out_mp}.')

        # Thermal boundary condition at input interface
        self.thermal_bc = None
        if "temperature" in self.input_variables:
            self.thermal_bc = "temperature"
        elif "heat_flux" in self.input_variables:
            self.thermal_bc = "heat_flux"

        # Phase change specific settings
        self.pc_settings = self.settings['PC']
        self.pcm_name = self.pc_settings.get('pcm_name', None)
        if self.pcm_name is None and self.multiphase:
            raise ValueError('PCM name in Fluent should be included in the json file in case of a multiphase simulation.')
        elif self.pcm_name is None and not self.multiphase:
            self.pcm_name = 'pcm'
        self.ini_condition = self.pc_settings.get('ini_condition', None)  # initial condition for outgoing variables
        self.melt_temp = self.pc_settings.get('melt_temp', None)
        if self.melt_temp is None:
            raise ValueError('Melt temperature should be provided in the json file in case of phase change problems')
        self.melt_enthalpy = self.pc_settings.get('melt_enthalpy', 0.0)
        if self.melt_enthalpy == 0.0:
            tools.print_info('No melting enthalpy is given: constant cp assumed for liquid solver.', layout='warning')
        self.volume_change = self.pc_settings.get('volume_change', True) # Account for volume change during melting
        if self.volume_change:
            tools.print_info('Density difference between solid and liquid accounted for during melting and rigid body motion.', layout='info')
        else:
            tools.print_info('Density difference between solid and liquid accounted for rigid body motion, not during melting.', layout='info')
        self.liquid_density = self.pc_settings.get('liquid_density', None)  # Give zero option and assume rho_s == rho_l
        if self.liquid_density is None:
            raise ValueError('Liquid density should be provided in the json file in case of unconstrained phase change problems.')
        self.solid_density = self.pc_settings.get('solid_density', 0.0) # Give zero option and assume rho_s == rho_l
        if self.solid_density == 0.0:
            self.solid_density == self.liquid_density.copy()
            tools.print_info('Equal density for solid and liquid is assumed.', layout='warning')

        # Rigid body motion specific settings
        self.rb_settings = self.settings['RB']
        self.restart_rb_only = self.rb_settings.get('restart_rb', 0)
        if self.restart_rb_only != 0 and self.restart:
            self.restart_rb_only = 0
            tools.print_info('Rigid body only restart overruled by general restart.', layout='warning')
        self.rb_tol = self.rb_settings.get('tolerance', 1E-6)
        self.rb_iter_min = self.rb_settings.get('iteration_min', 1)
        self.rb_iter_max = self.rb_settings.get('iteration_max', 10)
        if self.rb_iter_min > self.rb_iter_max:
            self.rb_iter_min = 1
            tools.print_info('rb_iter_min is set larger than rb_iter_max, rb_iter_min has been reset to 1.', layout='warning')
        self.rb_relax = self.rb_settings.get('relaxation', 1.0)  # default no relaxation
        self.fict_coeff = self.rb_settings.get('fict_coeff', 0.0)  # default no fictitious impedance
        self.multiplier = self.rb_settings.get('fict_multiplier', 1.0)  # default no multiplier for first timestep
        self.dyn_visc = self.rb_settings.get('liquid_dyn_visc', 0.0)  # Dynamic viscosity of liquid PCM
        self.fm_wall = self.rb_settings.get('fm_wall', None)
        if self.fict_coeff != 0.0 and self.dyn_visc == 0.0:
            raise ValueError('Liquid dyn. viscosity should be provided in the json file in case of fictitious damping.')
        self.rb_relax_method = self.rb_settings.get('relaxation_method', 'static') # or 'aitken'
        self.relax_first_ts =self.rb_settings.get('relax_first_ts', False)
        self.rot_update = self.rb_settings.get('rotational_update', 'quaternion') # 'off', 'rot_mat' or 'quaternion'
        self.rb_predictor = self.rb_settings.get('predictor', 'constant')  # or 'linear'
        self.buoyancy = self.rb_settings.get('buoyancy', True)
        self.weight_ramp = self.rb_settings.get('weight_ramp', 0)  # Nr. of time steps over which the full weight will be added
        self.gravity = self.rb_settings.get('gravity', [0, -9.81, 0])
        if self.restart_rb_only != 0:
            self.weight_ramp = 0
            tools.print_info('Weight ramp disabled due to rigid body only restart.', layout='warning')

        # contact model specific settings
        self.contact_settings = self.settings['contact_model']
        self.include_contact_force = self.contact_settings.get('contact_force', False)
        self.h_ul = self.contact_settings.get('upper_limit', 2e-4)  # [m] Upper limit
        self.h_ll = self.contact_settings.get('lower_limit', 2e-4)  # [m] Lower limit
        self.damping_ratio = self.contact_settings.get('damping_ratio', 1.0)
        self.k_mass = self.contact_settings.get('k_mass', 1.0)
        gap_walls = self.contact_settings.get('gap_walls', [])  # List of heated walls to calculate the gap from
        self.gap_ids = {}  # thread IDs corresponding to close contact walls
        self.gap_trees = {}  # thread IDs corresponding to close contact walls
        for gap_name in gap_walls:
            self.gap_ids[gap_name] = None
            self.gap_trees[gap_name] = None

        # Consistency checks
        if self.gap_ids:
            if self.fm_wall is None:
                raise ValueError('One wall from gap walls must be chosen as wall to calculate fictitious mass from.')
            elif self.fm_wall not in self.gap_ids:
                raise ValueError('Given wall to calculate fictitious mass from must be present within gap walls.')
        else:
            if self.fm_wall:
                raise ValueError('Wall given to calculate fictitious mass from but gap walls is empty.')
        if self.fict_coeff != 0.0 and not self.gap_ids:
            self.fict_coeff = 0.0
            tools.print_info('No gap wall given to calculate liquid layer thickness. Disabling fictitious impedance method.', layout='warning')

        # Initialise (with or without restart)
        if self.restart or self.restart_rb_only != 0:
            self.load_restart_rb_data()
        else:
            # Initialise rigid body motion kinematics
            self.v_trans_prev = np.zeros(3) # last time step translational velocity
            self.v_trans = np.zeros(3) # current translational velocity
            self.omega_prev = np.zeros(3) # last time step rotational velocity
            self.omega = np.zeros(3) # current rotational velocity
            self.a_trans_prev = np.zeros(3)  # last time step translational acceleration
            self.a_trans = np.zeros(3)  # current translational acceleration
            self.a_rot_prev = np.zeros(3)  # last time step rotational acceleration
            self.a_rot = np.zeros(3)  # current rotational acceleration
            self.com = np.zeros(3) # last time step center of mass
            self.com_prev = np.zeros(3) # current time step center of mass
            self.prev_patches = []
            self.new_patches = []
            self.force_pr_it = np.zeros(3)
            self.moment_pr_it = np.zeros(3)

        if not self.restart:
            if self.rot_update == 'rot_mat':
                self.orientation_prev = np.identity(3) # previous time step orientation matrix
                self.orientation = np.identity(3) # new time step orientation matrix
            else:
                self.orientation_prev = np.array([1.0, 0.0, 0.0, 0.0]) # previous time step orientation quaternion
                self.orientation = np.array([1.0, 0.0, 0.0, 0.0]) # new time step orientation quaternion

        # Initialise other variables
        self.volume = 0.0
        self.M_sys = 0.0
        self.a_trans_prev_it = np.zeros(3)
        self.force_int = np.zeros(3)
        self.moment_int = np.zeros(3)
        self.avg_v_trans = np.zeros(3)
        self.avg_omega = np.zeros(3)
        self.h_min = 0.0
        self.contact_patches = []

        # Aitken relaxation state variables
        if self.rb_relax_method == 'aitken':
            self.force_res_prev = np.zeros(3)
            self.moment_res_prev = np.zeros(3)
            self.aitken_relax_factor_force = self.rb_relax
            self.aitken_relax_factor_moment = self.rb_relax

        self.cnt = 0

    @tools.time_initialize
    def initialize(self):
        super().initialize()

        # prepare Fluent journal
        journal = f'v{self.version}.jou'
        thread_names_str = ''
        for thread_name in self.thread_ids:
            thread_names_str += ' "' + thread_name + '"'
        overset_thread_names_str = ''
        for thread_name in self.overset_thread_ids:
            overset_thread_names_str += ' "' + thread_name + '"'
        thread_domains_str = ''
        for thread_domain in self.settings['thread_domains']:
            thread_domains_str += ' "' + thread_domain + '"'
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
                    line = line.replace('|OVERSET_THREAD_NAMES|', overset_thread_names_str)
                    line = line.replace('|THREAD_DOMAINS|', thread_domains_str)
                    line = line.replace('|GAP_WALL|', '#t' if self.gap_ids else '#f')
                    line = line.replace('|UNSTEADY|', '#t' if self.unsteady else '#f')
                    line = line.replace('|MULTIPHASE|', '#t' if self.multiphase else '#f')
                    line = line.replace('|THERMAL_BC|', thermal_bc)
                    line = line.replace('|MATERIAL|', self.pcm_name)
                    line = line.replace('|FLOW_ITERATIONS|', str(self.flow_iterations))
                    line = line.replace('|DELTA_T|', str(self.delta_t))
                    line = line.replace('|TIMESTEP_START|', str(self.timestep_start))
                    line = line.replace('|END_OF_TIMESTEP_COMMANDS|', self.settings.get('end_of_timestep_commands','\n'))
                    line = line.replace('|END_OF_SETUP_COMMANDS|', self.pc_settings.get('end_of_setup_commands','\n'))
                    outfile.write(line)

        # prepare Fluent UDF
        if self.volume_change and self.solid_density != self.liquid_density:
            solid_density = self.solid_density
        else:
            solid_density = 0.0
        udf = 'udf_thermal.c'
        with open(join(self.dir_src, udf)) as infile:
            with open(join(self.dir_cfd, udf), 'w') as outfile:
                for line in infile:
                    line = line.replace('|MAX_NODES_PER_FACE|', str(self.mnpf))
                    line = line.replace('|TMP_DIRECTORY_NAME|', self.tmp_dir_unique)
                    line = line.replace('|TIME_STEP_SIZE|', str(self.delta_t))
                    line = line.replace('|MELT_TEMP|', str(self.melt_temp))
                    line = line.replace('|MELT_ENTHALPY|', str(self.melt_enthalpy))
                    line = line.replace('|SOLID_DENSITY|', str(solid_density))
                    line = line.replace('|LIQUID_DENSITY|', str(self.liquid_density))
                    line = line.replace('|TIME_STEP_START|', str(self.timestep_start))
                    line = line.replace('|UNSTEADY|', 'true' if self.unsteady else 'false')
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

        # pass on process to coco_messages for polling
        self.coco_messages.set_process(self.fluent_process)

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
        overset_names_found = []
        gap_walls_found = []
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
                    elif name in self.overset_thread_ids and name not in overset_names_found:
                        self.overset_thread_ids[name] = thread_id
                        overset_names_found.append(name)
                    elif name in self.gap_ids and name not in gap_walls_found:
                        self.gap_ids[name] = thread_id
                        gap_walls_found.append(name)
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
        with open(join(self.dir_cfd, 'bcs_overset.txt'), 'w') as file:
            file.write(f'{len(overset_names_found)}\n')
            for name, id in self.overset_thread_ids.items():
                file.write(f'{name} {id}\n')
        with open(join(self.dir_cfd, 'gap_wall.txt'), 'w') as file:
            file.write(f'{len(gap_walls_found)}\n')
            for name, id in self.gap_ids.items():
                file.write(f'{name} {id}\n')
        self.coco_messages.send_message('thread_ids_written_to_file')

        # remove "report.sum" because the batch options to overwrite report files and case files conflict in some
        # versions of Fluent (2023R1)
        os.unlink(report)

        # import node and face information if no restart
        if not self.restart:
            self.coco_messages.wait_message('nodes_and_faces_stored')

        # create Model used for coupling
        self.model = data_structure.Model()

        if not self.restart:
            # create Model for internal rigid body motion (in liquid solver only)
            self.model_rb = data_structure.Model()
            self.rb_model_settings = []

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
            # only nodal variables for rigid body interface
            if "nodes" in mp_name and not self.restart:
                self.model_rb.create_model_part(mp_name, x0, y0, z0, ids)
                self.rb_model_settings.append({"model_part": mp_name, "variables": ["prev_disp", "prev_melting_disp", "new_disp", "disp_step", "disp_step_melting"]})

        # create output ModelParts (nodes or faces)
        for j, item in enumerate(self.settings['interface_output']):
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
                data = np.loadtxt(file_name, skiprows=1, ndmin=2)
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

            # create initial conditions at output interface
            if self.ini_condition is not None:
                if "faces" in mp_name:
                    self.output_ini_cond[mp_name] = np.ones((data.shape[0], 1))*self.ini_condition

        # create interfaces
        self.interface_input = data_structure.Interface(self.settings['interface_input'], self.model)
        self.interface_output = data_structure.Interface(self.settings['interface_output'], self.model)
        if not self.restart:
            # create internal interface for rigid body motion
            self.interface_rb = data_structure.Interface(self.rb_model_settings, self.model_rb)

        # set initial conditions at output interface
        for key in self.output_ini_cond:
            if "faces" in key:
                for var in self.output_variables:
                    if var == 'heat_flux' or var == 'temperature':
                        if self.ini_condition is None:
                            raise ValueError('Initial condition must be defined in JSON for temperature or heat flux')
                        else:
                            self.interface_output.set_variable_data(key, var, self.output_ini_cond[key])

        # Check rigid body interface in case of restart
        if self.restart:
            if not self.interface_input.has_same_model_parts(self.interface_rb):
                raise ValueError('Restart not possible because model parts in reloaded rigid body interface do not match those of the new input interface.')

        if self.gap_ids:
            for wall_name in self.gap_ids:
                try:
                    gap_file = join(self.dir_cfd, f'nodes_{wall_name}.dat')
                    df = pd.read_csv(gap_file, delim_whitespace=True, skiprows=1, header=None)

                    # Extract coords (2D or 3D) and remove duplicates
                    raw = df.iloc[:, 0:2].values if df.shape[1] == 3 else df.iloc[:, 0:3].values
                    coords = np.unique(raw, axis=0)

                    wall_tree = cKDTree(coords)
                    print(f"KD-Tree initialized with {len(coords)} nodes for wall {wall_name}.")
                    self.gap_trees[wall_name] = wall_tree

                    # --- Calculate Length of Correct Wall (Greedy Walk) ---
                    if wall_name == self.fm_wall:
                        curr_idx = np.argmin(coords[:, 0])  # Start at min X
                        visited = {curr_idx}
                        l_current = 0.0

                        while len(visited) < len(coords):
                            # Query neighbors (k=10 buffer for visited nodes)
                            dists, idxs = wall_tree.query(coords[curr_idx], k=10)

                            # Find nearest unvisited neighbor
                            next_step = next(((d, i) for d, i in zip(dists, idxs) if i not in visited), None)

                            if next_step:
                                dist, idx = next_step
                                l_current += dist
                                visited.add(idx)
                                curr_idx = idx
                            else:
                                print(f"Warning: Discontinuity in wall {wall_name}")
                                break

                        self.wall_length = l_current
                        print(f"Wall Length: {self.wall_length:.6f} [m]")

                except Exception as e:
                    print(f"Error loading wall nodes: {e}")

    def initialize_solution_step(self):
        super().initialize_solution_step()

        self.iteration = 0
        self.timestep += 1

        # save for linear predictor
        if self.rb_predictor == 'linear':
            v_trans_prev2 = self.v_trans_prev.copy()
            omega_prev2 = self.omega_prev.copy()

        # update rigid body kinematics
        self.v_trans_prev = self.v_trans.copy()
        self.omega_prev = self.omega.copy()
        self.a_trans_prev = self.a_trans.copy()
        self.a_rot_prev = self.a_rot.copy()
        self.orientation_prev = self.orientation.copy()
        self.com_prev = self.com.copy()
        self.prev_patches = self.new_patches.copy()

        # Reset Aitken relaxation residuals for the new time step
        if self.rb_relax_method == 'aitken':
            self.force_res_prev = np.zeros(3)
            self.moment_res_prev = np.zeros(3)
            self.aitken_relax_factor_force = np.clip(self.aitken_relax_factor_force, 0.01, self.rb_relax)
            self.aitken_relax_factor_moment = np.clip(self.aitken_relax_factor_moment, 0.01, self.rb_relax)

        # linear predictor step for v_trans & omega_z
        if self.rb_predictor == 'linear' and self.timestep > 1:
            self.v_trans = 2 * self.v_trans_prev - v_trans_prev2
            self.omega = 2 * self.omega_prev - omega_prev2
            # Make accelerations consistent
            self.a_trans = 2 * (self.v_trans - self.v_trans_prev) / self.delta_t - self.a_trans_prev
            self.a_rot = 2 * (self.omega - self.omega_prev) / self.delta_t - self.a_rot_prev

        # Update previous displacement
        for item in self.settings['interface_input']:
            mp_name = item['model_part']
            if "nodes" in mp_name:
                self.interface_rb.set_variable_data(mp_name, "prev_disp", self.interface_rb.get_variable_data(mp_name, "new_disp"))
                self.interface_rb.set_variable_data(mp_name, "prev_melting_disp", self.interface_input.get_variable_data(mp_name, "displacement"))

        self.coco_messages.send_message('next')
        self.coco_messages.wait_message('next_ready')

    @tools.time_solve_solution_step
    def solve_solution_step(self, interface_input):
        self.rb_iter = 0
        check_tolerance = False # rigid body iteration convergence
        self.iteration += 1

        # process input interface data
        # store incoming variables
        self.interface_input.set_interface_data(interface_input.get_interface_data()) # total melting displacement

        for var in self.input_variables:
            if var == "displacement":
                # calculate single time step displacement due to melting
                self.melting_displacement()
            else:
                # write interface input data
                self.write_input_to_file(var)

                # copy input data for debugging
                if self.debug:
                    for dct in self.interface_input.parameters:
                        mp_name = dct['model_part']
                        thread_id = self.model_part_thread_ids[mp_name]
                        src = accepted_variables_pc_liquid_rb['in'][var][0] + f'_timestep{self.timestep}_thread{thread_id}.dat'
                        dst = accepted_variables_pc_liquid_rb['in'][var][0] + f'_timestep{self.timestep}_thread{thread_id}_Iter{self.iteration}.dat'
                        cmd = f'cp {join(self.dir_cfd, src)} {join(self.dir_cfd, dst)}'
                        os.system(cmd)

        # Start new coupling iteration
        self.coco_messages.send_message('continue')
        self.coco_messages.wait_message('continue_ready')

        self.print_rb_header()

        # Rigid body motion iterations
        while (not check_tolerance or self.rb_iter < self.rb_iter_min) and self.rb_iter < self.rb_iter_max:
            self.rb_iter += 1

            self.update_nodal_positions()
            self.update_file_rb_motion()  # necessary for momentum source term and overset boundary motion

            force_last_it = self.force_int
            moment_last_it = self.moment_int

            # let Fluent run
            self.coco_messages.send_message('rigidbody')
            self.coco_messages.wait_message('rigidbody_ready')

            if self.debug:
                src = f'rigid_body_timestep{self.timestep}.dat'
                dst = f'rigid_body_timestep{self.timestep}_Iter{self.iteration}.dat'
                cmd = f'cp {join(self.dir_cfd, src)} {join(self.dir_cfd, dst)}'
                os.system(cmd)

            # Read new rigid body displacement
            self.rigid_body_motion()

            # check convergence of translational and rotational velocity
            res_f = np.linalg.norm(self.force_int - force_last_it)
            res_m = np.linalg.norm(self.moment_int - moment_last_it)
            check_tolerance = (res_f <= self.rb_tol) and (res_m <= self.rb_tol)

            # print rigid body iteration information to terminal
            self.print_rb_iteration_info(res_f, res_m)

            # only one iteration for the 1st coupling iteration --> avoids divergence
            if not self.relax_first_ts and self.timestep == 1 and self.iteration == 1:
                break

        self.cnt += self.rb_iter

        # print warning when tolerance not reached within maximum number of iterations
        if self.rb_iter_max == 1:
            tools.print_info('Explicit rigid body update: no iterations until convergence.')
        if self.rb_iter >= self.rb_iter_max and self.rb_iter_max > 1:
            if np.linalg.norm(self.force_int - force_last_it) > self.rb_tol:
                warning_text = f'Rigid body iterations did not converge below tolerance for force vector: {np.linalg.norm(self.force_int - force_last_it)} > {self.rb_tol}.'
                tools.print_info(warning_text, layout='warning')
            if np.linalg.norm(self.moment_int - moment_last_it) > self.rb_tol:
                warning_text = f'Rigid body iterations did not converge below tolerance for moment: {np.linalg.norm(self.moment_int - moment_last_it)} > {self.rb_tol}.'
                tools.print_info(warning_text, layout='warning')

        # process output interface data once rigid body iterations are converged
        for dct in self.interface_output.parameters:
            mp_name = dct['model_part']
            thread_id = self.model_part_thread_ids[mp_name]

            # read in datafile
            for var in dct['variables']:
                prefix = accepted_variables_pc_liquid_rb['out'][var][0]
                # Avoid repeat of commands in case variables are stored in the same file

                if var == 'displacement':
                    if 'nodes' not in mp_name:
                        raise ValueError('Model part must be node-based for the displacement variable')

                    # Pass input displacement directly to output for use in solid solver
                    self.interface_output.set_variable_data(mp_name, 'displacement',
                                                            self.interface_input.get_variable_data(mp_name.replace('out', 'in'), 'displacement'))

                elif prefix != 'skip':
                    data = self.read_output_file(prefix, thread_id)
                    if accepted_variables_pc_liquid_rb['out'][var][1] == 0:
                        req_dim = self.dimensions + 1 + self.mnpf
                    else:
                        req_dim = min([self.dimensions, accepted_variables_pc_liquid_rb['out'][var][1]]) + self.mnpf
                    if data.shape[1] != req_dim:
                        raise ValueError('Given dimension does not match coordinates')

                    if accepted_variables_pc_liquid_rb['out'][var][1] == 0:
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

                    if accepted_variables_pc_liquid_rb['out'][var][1] == 3:
                        # get face coordinates and ids
                        vector_tmp = np.zeros((data.shape[0], 3)) * 0.
                        vector_tmp[:, :self.dimensions] = data[:, : -self.mnpf]
                        ids_tmp, _ = self.get_unique_face_ids(data[:, -self.mnpf:])

                        # sort and remove doubles
                        args = np.unique(ids_tmp, return_index=True)[1].tolist()
                        vector = vector_tmp[args, :]
                        ids = ids_tmp[args]

                        # store vector values in interface
                        model_part = self.model.get_model_part(mp_name)
                        if ids.size != model_part.size:
                            raise ValueError('Size of data does not match size of ModelPart')
                        if not np.all(ids == model_part.id):
                            raise ValueError('IDs of data do not match ModelPart IDs')

                        self.interface_output.set_variable_data(mp_name, var, vector)

                    if accepted_variables_pc_liquid_rb['out'][var][1] == 1:
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

        # update a report-file each timestep with the RB values
        self.update_report_file()

    @tools.time_save
    def output_solution_step(self):
        super().output_solution_step()

        # save data for restart
        if self.save_restart != 0 and self.timestep % self.save_restart == 0:
            self.save_restart_rb_data()
            if self.save_restart < 0 and self.timestep + self.save_restart > self.timestep_start_current:
                try:
                    os.remove(self.case_name + f'_restart_ts{self.timestep + self.save_restart}.pickle')
                except OSError:
                    pass

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
        print(f"Average nr. of RB iterations = {self.cnt/self.timestep}")
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
                try:
                    os.remove(join(self.dir_cfd, f'nodes_pc_timestep{timestep}_thread{thread_id}.dat'))
                except OSError:
                    pass
                try:
                    os.remove(join(self.dir_cfd, f'rigid_body_timestep{timestep}.dat'))
                except OSError:
                    pass
                for var in self.output_variables:
                    if accepted_variables_pc_liquid_rb['out'][var][0] != 'skip':
                        try:
                            os.remove(join(self.dir_cfd, accepted_variables_pc_liquid_rb['out'][var][0] + f'_timestep{timestep}_thread{thread_id}.dat'))
                        except OSError:
                            pass
                for var in self.input_variables:
                    try:
                        os.remove(join(self.dir_cfd, accepted_variables_pc_liquid_rb['in'][var][0] + f'_timestep{timestep}_thread{thread_id}.dat'))
                    except OSError:
                        pass
            try:
                os.remove(join(self.dir_cfd, f'RB_update_timestep{self.timestep}.dat'))
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
        if var == 'temperature' or var == 'heat_flux':
            for dct in self.interface_input.parameters:
                mp_name = dct['model_part']
                if 'faces' in mp_name:
                    thread_id = self.model_part_thread_ids[mp_name]
                    model_part = self.model.get_model_part(mp_name)
                    scalar = self.interface_input.get_variable_data(mp_name, var)
                    face_nodeIDs = np.zeros((np.size(model_part.id), self.mnpf))
                    for index, id in enumerate(model_part.id):
                        list = self.reverse_face_ids(id, mp_name)
                        if len(list) == self.mnpf:
                            face_nodeIDs[index,:] = np.array(list)
                        else:
                            raise ValueError("Nr. of nodes returned from reverse_face_ids function does not correspond to mnpf")
                    prof = np.append(scalar, face_nodeIDs, axis=1)
                    fmt = '%27.17e'
                    for i in range(self.mnpf):
                        fmt += '%27d'
                    prefix = accepted_variables_pc_liquid_rb['in'][var][0]
                    tmp = prefix + f'_timestep{self.timestep}_thread{thread_id}.dat'
                    file_name = join(self.dir_cfd, tmp)
                    np.savetxt(file_name, prof, fmt=fmt, header='temperature unique-ids', comments='')

    def update_nodal_positions(self):
        if len(self.interface_input.parameters) > 1:
            if self.timestep == 1 and self.iteration == 1 and self.rb_iter == 1:
                tools.print_info('Solver wrapper currently only supports one rigid body. Only one model part will receive nodal updates by update_nodal_positions().', layout='warning')
        for j, dct in enumerate(self.interface_input.parameters):
            if j == 0:
                mp_name = dct['model_part']
                if 'nodes' in mp_name:
                    thread_id = self.model_part_thread_ids[mp_name]
                    model_part = self.model.get_model_part(mp_name)
                    last_full_disp = self.interface_rb.get_variable_data(mp_name, 'prev_disp')
                    disp_step_melting = self.interface_rb.get_variable_data(mp_name, 'disp_step_melting')
                    x_prev = model_part.x0 + last_full_disp[:, 0]
                    y_prev = model_part.y0 + last_full_disp[:, 1]
                    z_prev = model_part.z0 + last_full_disp[:, 2]
                    r_prev = np.column_stack((x_prev, y_prev, z_prev))

                    # update r_prev with the melting displacement, rotated to the global frame (liquid domain)
                    r_temp = r_prev + disp_step_melting

                    # Crank-Nicolson (trapezoidal) integration
                    # Use the average velocity over the time step for a second-order accurate position update
                    self.avg_v_trans = 0.5 * (self.v_trans_prev + self.v_trans)
                    self.avg_omega = 0.5 * (self.omega_prev + self.omega)

                    omega_array = np.zeros(np.shape(r_prev))
                    omega_array[:, 2] = self.avg_omega[2]
                    disp_step_rb = (self.avg_v_trans + np.cross(omega_array, (r_temp - self.com_prev))) * self.delta_t

                    self.interface_rb.set_variable_data(mp_name, 'new_disp', last_full_disp + disp_step_melting + disp_step_rb)

                    r_new = r_temp + disp_step_rb
                    x = r_new[:, 0]
                    y = r_new[:, 1]
                    z = r_new[:, 2]

                    if self.dimensions == 2:
                        data = np.rec.fromarrays([x, y, model_part.id])
                        fmt = '%27.17e%27.17e%27d'
                    else:
                        data = np.rec.fromarrays([x, y, z, model_part.id])
                        fmt = '%27.17e%27.17e%27.17e%27d'
                    tmp = f'nodes_update_timestep{self.timestep}_thread{thread_id}.dat'
                    file_name = join(self.dir_cfd, tmp)
                    np.savetxt(file_name, data, fmt=fmt, header=f'{model_part.size}', comments='')

    def update_file_rb_motion(self):
        if len(self.interface_input.parameters) > 1:
            if self.timestep == 1 and self.iteration == 1 and self.rb_iter == 1:
                tools.print_info('Solver wrapper currently only supports one rigid body. Only one "RB_update" file will be written by update_file_rb_motion().', layout='warning')
        for j, dct in enumerate(self.interface_input.parameters):
            if j == 0:
                mp_name = dct['model_part']
                overset_boundary_name = None
                if 'nodes' in mp_name:
                    for key, value in self.outer_inner_surf_dict.items():
                        if key in mp_name:  # check substring match
                            overset_boundary_name = value
                            break  # stop at the first match

                    if overset_boundary_name is None and self.define_os:
                        raise ValueError(f"No overset boundary mapping found for model part '{mp_name}'")

                    # thread_id = self.overset_thread_ids[overset_boundary_name]
                    omega = np.zeros(3)
                    omega[2] = self.avg_omega[2]

                    tmp = f'RB_update_timestep{self.timestep}.dat'
                    file_name = join(self.dir_cfd, tmp)

                    with open(file_name, "w") as f:
                        f.write(" ".join(f"{d:.16e}" for d in self.avg_v_trans) + "\n")
                        f.write(" ".join(f"{d:.16e}" for d in omega) + "\n")
                        f.write(" ".join(f"{d:.16e}" for d in self.com_prev) + "\n")

    def melting_displacement(self):
        for dct in self.interface_input.parameters:
            mp_name = dct['model_part']
            if 'nodes' in mp_name:
                thread_id = self.model_part_thread_ids[mp_name]
                model_part = self.model.get_model_part(mp_name)
                total_melting_disp = self.interface_input.get_variable_data(mp_name, 'displacement')
                prev_melting_disp = self.interface_rb.get_variable_data(mp_name, 'prev_melting_disp')

                # Single time step displacement due to melting (in the frame of reference of the solid)
                disp_step_melting = total_melting_disp - prev_melting_disp

                # Write the melting only coordinates to a file in order to calculate the local volume change in the mesh
                x_melting = model_part.x0 + total_melting_disp[:, 0]
                y_melting = model_part.y0 + total_melting_disp[:, 1]
                z_melting = model_part.z0 + total_melting_disp[:, 2]

                if self.dimensions == 2:
                    data = np.rec.fromarrays([x_melting, y_melting, model_part.id])
                    fmt = '%27.17e%27.17e%27d'
                else:
                    data = np.rec.fromarrays([x_melting, y_melting, z_melting, model_part.id])
                    fmt = '%27.17e%27.17e%27.17e%27d'
                tmp = f'nodes_pc_timestep{self.timestep}_thread{thread_id}.dat'
                file_name = join(self.dir_cfd, tmp)
                np.savetxt(file_name, data, fmt=fmt, header=f'{model_part.size}', comments='')

                if self.rot_update == 'rot_mat':
                    # Convert from the local frame (solid) to the global frame of the rigid body in the liquid domain
                    disp_step_melting_global = np.dot(disp_step_melting, self.orientation_prev.T)
                else:
                    # Use the quaternion from the start of the step to rotate the displacement vectors
                    disp_step_melting_global = quat_rotate_vector_array(self.orientation_prev, disp_step_melting)

                self.interface_rb.set_variable_data(mp_name, 'disp_step_melting', disp_step_melting_global)

    def rigid_body_motion(self):
        tmp = f'rigid_body_timestep{self.timestep}.dat'
        file_name = join(self.dir_cfd, tmp)
        with open(file_name, 'r') as f:
            lines = [line.strip() for line in f if line.strip()]  # remove blank lines

        # Skip the header line
        self.volume = float(lines[1])
        self.com = np.array(lines[2].split(), dtype=float)
        if self.timestep == 1 and self.iteration == 1:
            self.com_prev = self.com

        moi = np.array(lines[3].split(), dtype=float)
        # fluid-integrated forces read from CFD (these are the ones we will relax)
        self.force_int = np.array(lines[4].split(), dtype=float)
        self.moment_int = np.array(lines[5].split(), dtype=float)

        # Keep a copy of the raw fluid force and moment read from file (before contact)
        force_raw = self.force_int.copy()
        moment_raw = self.moment_int.copy()

        # --- FORCE RELAXATION (static or aitken) ---
        # Option to skip relaxation on the very first timestep and iteration
        if not self.relax_first_ts and self.timestep == 1 and self.iteration == 1:
            force_relaxed = force_raw.copy()
        else:
            if self.rb_relax_method == 'static':
                # simple under-relaxation on forces
                force_relaxed = self.force_pr_it + self.rb_relax * (force_raw - self.force_pr_it)
            elif self.rb_relax_method == 'aitken':
                # Aitken on the force residuals
                self.aitken_relax_factor_force = self.calc_aitken(
                    force_raw, self.force_pr_it, self.force_res_prev, self.aitken_relax_factor_force
                )
                # update previous residual for next aitken step
                self.force_res_prev = force_raw - self.force_pr_it
                force_relaxed = self.force_pr_it + self.aitken_relax_factor_force * (force_raw - self.force_pr_it)
            else:
                # fallback: no relaxation
                force_relaxed = force_raw.copy()

        # store relaxed force as "previous" for next iteration
        self.force_pr_it = force_relaxed.copy()

        # --- MOMENT RELAXATION (static or aitken) ---
        # Option to skip relaxation on the very first timestep and iteration
        if not self.relax_first_ts and self.timestep == 1 and self.iteration == 1:
            moment_relaxed = moment_raw.copy()
        else:
            if self.rb_relax_method == 'static':
                # simple under-relaxation on moments
                moment_relaxed = self.moment_pr_it + self.rb_relax * (moment_raw - self.moment_pr_it)
            elif self.rb_relax_method == 'aitken':
                # Aitken on the moment residuals
                self.aitken_relax_factor_moment = self.calc_aitken(
                    moment_raw, self.moment_pr_it, self.moment_res_prev, self.aitken_relax_factor_moment
                )
                # update previous residual for next aitken step
                self.moment_res_prev = moment_raw - self.moment_pr_it
                moment_relaxed = self.moment_pr_it + self.aitken_relax_factor_moment * (moment_raw - self.moment_pr_it)
            else:
                # fallback: no relaxation
                moment_relaxed = moment_raw.copy()

        # store relaxed moment as "previous" for next iteration
        self.moment_pr_it = moment_relaxed.copy()

        # Calculate fraction of weight that is accounted for
        self.weight_factor = self.timestep / self.weight_ramp if self.timestep < self.weight_ramp else 1

        # --- TRANSLATIONAL UPDATE ---
        g = np.array(self.gravity)  # gravitational acceleration
        mass_solid = self.solid_density * self.volume

        # Fictitious Mass / Damping method to anticipate high added mass & viscous damping in fluid solver
        if self.gap_ids:
            h_used, self.h_min = self.calculate_h_min()

            print("\n")
            print(f'Min. gap width = {self.h_min * 1000} mm')
            print(f'Number of active contact patches = {len(self.contact_patches)}')

            # Fictitious mass components are always calculated at the chosen wall
            epsilon = 1e-4 * self.wall_length
            h_eff = max(h_used, epsilon)  # Singularity protection
            depth = 1.0 # in 2D

            fict_mass = self.liquid_density * (self.wall_length ** 3) * depth / h_eff
            fict_damping = depth * (self.wall_length ** 3) * (self.solid_density / self.liquid_density) * self.dyn_visc / (h_eff ** 3)
        else:
            fict_mass = 0
            fict_damping = 0

        self.a_trans_prev_it = self.a_trans.copy()

        if self.buoyancy:
            F_net = force_relaxed + self.weight_factor * (self.solid_density - self.liquid_density) * g * self.volume
        else:
            F_net = force_relaxed + self.weight_factor * mass_solid * g

        # --- PREPARE CONTACT FORCES ---
        contact_F = np.zeros(3)
        contact_M = np.zeros(3)

        if self.include_contact_force and self.gap_ids:
            contact_F, contact_M = self.calc_contact_force_and_moment()

        # Add to Net Force
        F_net += contact_F

        # We add (M_sys * a_prev_iter) to the forces
        # This dampens the acceleration update significantly without altering the final converged physics.
        self.M_sys = fict_mass + self.delta_t * fict_damping / 2
        M_sys_raw = self.M_sys.copy()
        self.M_sys *= self.fict_coeff
        if self.rb_iter == 1:
            self.M_sys *= self.multiplier
        F_tot = F_net + self.M_sys * self.a_trans_prev_it

        # Update velocity with Crank-Nicolson integration (2nd order)
        self.a_trans = F_tot / (mass_solid + self.M_sys)
        self.v_trans = self.v_trans_prev +  0.5 * (self.a_trans_prev + self.a_trans) * self.delta_t

        # --- TRANSLATIONAL WALL GUARD (FAILSAFE) ---
        if self.h_min < self.h_ll and self.contact_patches:
            deepest_contact = self.contact_patches[0][0]
            normal = deepest_contact['normal']
            if np.dot(self.v_trans, normal) < 0:
                print(f"--- HARD STOP ACTIVATED ---")
                rejected = np.dot(self.v_trans, normal) * normal
                self.v_trans = self.v_trans - rejected
                self.a_trans = 2 * (self.v_trans - self.v_trans_prev) / self.delta_t - self.a_trans_prev

        # --- ROTATIONAL UPDATE ---
        if self.rot_update == 'off':
            # Enforce zero rotation when turned off
            self.a_rot = np.zeros(3)
            self.omega = np.zeros(3)
            self.orientation = np.array([1.0, 0.0, 0.0, 0.0]) # Identity quaternion
        else:
            # NOTE: This implementation assumes a DIAGONAL moment of inertia tensor,
            # where the 'moi' vector represents [I_xx, I_yy, I_zz]. For a full
            # 3x3 tensor, a matrix inversion would be required.
            # Calculate raw rotational acceleration
            fict_moi = (self.M_sys / (self.solid_density * self.volume)) * moi
            moi_total = moi + fict_moi
            a_rot_prev_it = self.a_rot.copy()
            # moi is in this case a vector representing the diagonal elements of the otherwise empty moi matrix
            moment_tot = moment_relaxed + contact_M + fict_moi * a_rot_prev_it

            self.a_rot = np.divide(moment_tot, moi_total, out=np.zeros_like(moment_tot), where=moi_total != 0)

            # Update rotational velocity with Crank-Nicolson integration (2nd order)
            self.omega = self.omega_prev + 0.5 * (self.a_rot_prev + self.a_rot) * self.delta_t

            if self.h_min < self.h_ll and self.contact_patches:
                deepest_contact = self.contact_patches[0][0]
                normal = deepest_contact['normal']
                pinch_point = deepest_contact['point']
                rot_vel = np.cross(self.omega, (pinch_point - self.com))
                if np.dot(rot_vel, normal) < 0:
                    self.omega = np.zeros_like(self.omega)
                    # Make acceleration consistent
                    self.a_rot = 2 * (self.omega - self.omega_prev) / self.delta_t - self.a_rot_prev

            if self.rot_update == 'quaternion':
                # Update orientation using quaternions for second-order accuracy
                avg_omega = 0.5 * (self.omega_prev + self.omega)
                delta_quat = quat_from_angular_velocity(avg_omega, self.delta_t)
                self.orientation = quat_multiply(delta_quat, self.orientation_prev)
                self.orientation = quat_normalize(self.orientation)
            elif self.rot_update == 'rot_mat':
                # Update rotation matrix
                avg_omega = 0.5 * (self.omega_prev + self.omega)
                theta = self.delta_t * avg_omega
                self.orientation = self.orientation_prev @ exp_map(theta)  # Only rotation along z-axis assumed

                # Re-orthonormalize if drift exceeds tolerance
                if np.linalg.norm(self.orientation.T @ self.orientation - np.identity(3)) > 1e-8:
                    U, _, Vt = np.linalg.svd(self.orientation)
                    self.orientation = U @ Vt
    
    def read_output_file(self, prefix, thread_id=None):
        tmp = prefix + f'_timestep{self.timestep}_thread{thread_id}.dat'
        file_name = join(self.dir_cfd, tmp)
        data = np.loadtxt(file_name, skiprows=1, ndmin=2)
        # copy output data for debugging
        if self.debug:
            dst = prefix + f'_timestep{self.timestep}_thread{thread_id}_it{self.iteration}.dat'
            cmd = f'cp {file_name} {join(self.dir_cfd, dst)}'
            os.system(cmd)
        return data

    def calculate_h_min(self):
        """
        Groups contact nodes into spatially distinct 'Patches'.
        Returns a list of patches, where each patch is a list of candidate dictionaries.
        """
        r = None
        # --- 1. Extract Interface Node Coordinates (Unchanged) ---
        for dct in self.interface_input.parameters:
            mp_name = dct['model_part']
            if 'nodes' in mp_name:
                model_part = self.model.get_model_part(mp_name)
                # Work with new displacement to update each coupling iteration!
                last_full_disp = self.interface_rb.get_variable_data(mp_name, 'new_disp')
                x = model_part.x0 + last_full_disp[:, 0]
                y = model_part.y0 + last_full_disp[:, 1]
                if self.dimensions == 3:
                    z = model_part.z0 + last_full_disp[:, 2]
                    r = np.column_stack((x, y, z))
                else:
                    r = np.column_stack((x, y))
                break

        if r is None:
            return None, None

        # --- 2. KDTree Query & Danger Zone Detection ---
        h_used = None
        h_min = None

        # We will collect ALL candidates from ALL walls first
        # Format: (distance, normal_vector, contact_point_coords)
        candidates = []

        for wall_name in self.gap_trees:
            wall_tree = self.gap_trees[wall_name]

            # Find nearest distance from any interface node to the wall
            dists, wall_ids = wall_tree.query(r, k=1, workers=-1)

            h = np.min(dists)
            h_min = min(h_min, h) if h_min is not None else h

            # Save h_used for Fictitious Mass (specific to fm_wall)
            if wall_name == self.fm_wall:
                h_used = h

            # Identify nodes inside the Danger Zone
            danger_mask = dists < self.h_ul
            danger_indices = np.where(danger_mask)[0]

            if len(danger_indices) > 0:
                close_dists = dists[danger_indices]
                close_ids_wall = wall_ids[danger_indices]
                close_p_itf = r[danger_indices]
                close_p_wall = wall_tree.data[close_ids_wall]

                # Calculate normals for these points
                diff_vecs = close_p_itf - close_p_wall
                norms = np.linalg.norm(diff_vecs, axis=1)
                norms[norms < 1e-12] = 1.0
                close_normals = diff_vecs / norms[:, None]

                # Handle 2D -> 3D conversion for normals/points if needed
                if self.dimensions == 2:
                    z_col = np.zeros((len(danger_indices), 1))
                    close_normals = np.hstack((close_normals, z_col))
                    close_p_itf = np.hstack((close_p_itf, z_col))

                for j in range(len(danger_indices)):
                    candidates.append({
                        'h': close_dists[j],
                        'normal': close_normals[j],
                        'point': close_p_itf[j],
                        'id': danger_indices[j]  # nice for debugging
                    })

        # --- 3. CLUSTERING (True Chaining / Flood Fill) ---
        candidates.sort(key=lambda x: x['h'])  # Deepest first

        self.contact_patches = []
        processed_indices = set()

        sep_tol = 2 * self.h_ul

        for cand in candidates:
            # If this node is already part of a chain, skip it
            if cand['id'] in processed_indices:
                continue

            # Start a NEW patch with this deepest node
            current_patch = [cand]
            processed_indices.add(cand['id'])

            # Initialize the search queue with the leader
            search_queue = [cand]

            # --- FLOOD FILL LOOP ---
            # Keep searching until we run out of connected neighbors
            while len(search_queue) > 0:
                # Pop the next node to expand from
                expansion_node = search_queue.pop(0)
                expansion_point = expansion_node['point']

                # Check ALL candidates to see if they are neighbors of 'expansion_node'
                for potential_neighbor in candidates:
                    # Skip if already processed
                    if potential_neighbor['id'] in processed_indices:
                        continue

                    # Calculate distance to the CURRENT expansion node (not just the leader)
                    dist = np.linalg.norm(potential_neighbor['point'] - expansion_point)

                    # If connected, add to patch AND to queue (to extend the chain further)
                    if dist < sep_tol:
                        processed_indices.add(potential_neighbor['id'])
                        current_patch.append(potential_neighbor)
                        search_queue.append(potential_neighbor)

            # The chain is exhausted, save the patch
            self.contact_patches.append(current_patch)

        if len(self.contact_patches) > 0:
            print(f"CONTACT: Found {len(self.contact_patches)} patch(es).")
            for idx, patch in enumerate(self.contact_patches):
                print(f"  Patch {idx}: {len(patch)} nodes")

        return h_used, h_min

    def calc_contact_force_and_moment(self):
        """
        Calculates repulsive force using Normalized Weighted Average (NWA)
        to smooth transitions between nodes.
        """
        total_force = np.zeros(3)
        total_moment = np.zeros(3)

        # 1. Define Stiffness (k) and Damping (c)
        # We want the "collision" to be resolved over roughly 10 timesteps.

        # Based on harmonic oscillator period T = 2*pi*sqrt(m/k)
        # We want a half-period (impact) to match impact_duration.
        # k = m * (pi / impact_duration)^2
        mass = self.k_mass * self.solid_density * self.volume
        resolution_steps = 10
        k_stiff = mass * (np.pi / (resolution_steps * self.delta_t)) ** 2

        # c_crit = 2 * sqrt(m * k)
        c_damp = self.damping_ratio * 2 * np.sqrt(mass * k_stiff)

        # 2. Prepare storage for NEXT step
        # We store: {'p_eff': vector, 'h_eff': float}
        available_prev = list(self.prev_patches)
        self.new_patches = []

        # Define a spatial tolerance to recognize "the same patch"
        match_tolerance = 2.0 * self.h_ul

        # 3. Loop over Patches
        for patch in self.contact_patches:

            # --- A. Accumulate Weighted Averages ---
            w_sum = 0.0
            weighted_normal = np.zeros(3)
            weighted_point = np.zeros(3)
            weighted_h = 0.0

            for node in patch:
                h_i = node['h']

                # Weight Function: Linear kernel
                w_i = max(0.0, self.h_ul - h_i) / self.h_ul

                w_sum += w_i
                weighted_normal += w_i * node['normal']
                weighted_point += w_i * node['point']
                weighted_h += w_i * h_i

            if w_sum <= 1e-12:
                continue

            # --- B. Normalize ---
            # 1. Average Normal (Direction)
            n_eff = weighted_normal / w_sum
            n_eff = n_eff / np.linalg.norm(n_eff)  # Re-normalize to unit vector

            # 2. Average Position (Center of Pressure)
            p_eff = weighted_point / w_sum

            # 3. Average Gap (Penetration Depth)
            h_eff = weighted_h / w_sum

            # --- C. TRACKING: Find the closest previous patch ---
            v_eff = 0.0

            # Check history if it exists
            if available_prev:
                best_dist = float('inf')
                best_index = -1

                # Search for the spatially closest patch in the available pool
                for i, prev in enumerate(available_prev):
                    dist = np.linalg.norm(p_eff - prev['p_eff'])
                    if dist < best_dist:
                        best_dist = dist
                        best_index = i

                # If the closest patch is valid, pop it
                if best_dist < match_tolerance and best_index != -1:
                    matched_patch = available_prev.pop(best_index)

                    # Calculate velocity
                    v_eff = (h_eff - matched_patch['h_eff']) / self.delta_t

            # Save current data for next step
            self.new_patches.append({'p_eff': p_eff, 'h_eff': h_eff})

            # --- D. Calculate Force on this Effective Contact ---
            # 1. Spring Force: Linear (n=1) or Hertzian (n=1.5)
            penetration = max(0.0, self.h_ul - h_eff)
            f_spring_mag = k_stiff * penetration
            print(f'Spring force = {f_spring_mag} N')

            # 2. Hunt-Crossley Damping
            f_damp_mag = -c_damp * 2 * (penetration / (self.h_ul - self.h_ll)) * v_eff
            print(f'v_eff = {v_eff * 1e6} µm/s')
            print(f'Damper force = {f_damp_mag} N')

            # 3. Total Patch Force
            f_mag = max(0.0, f_spring_mag + f_damp_mag)
            f_vec = f_mag * n_eff

            # 4. Add to Body Sums
            total_force += f_vec
            r_vec = p_eff - self.com
            total_moment += np.cross(r_vec, f_vec)

        return total_force, total_moment

    def update_report_file(self):
        """Create or update the rigid-body report file each timestep."""

        tmp = "rigid-body-report-file.out"
        file_name = join(self.dir_cfd, tmp)

        # If first timestep: create file and write header
        if self.timestep == 1:
            with open(file_name, "w") as f:
                f.write("# CoCoNuT rigid body motion history\n")
                f.write("#\n")
                f.write("#  {:>10}  {:>12}  {:>12}  {:>12}  {:>12}  {:>12}  {:>12}  {:>12}  {:>12}  {:>12}  {:>12}\n"
                        .format("time", "CG_X", "CG_Y", "V_X", "V_Y", "THETA_Z", "F_X", "F_Y", "M_Z", "volume", "h_min"))
                f.write("#  {:>10}  {:>12}  {:>12}  {:>12}  {:>12}  {:>12}  {:>12}  {:>12}  {:>12}  {:>12}  {:>12}\n"
                        .format("(s)", "(m)", "(m)", "(m/s)", "(m/s)", "(deg)", "(N)", "(N)", "(N*m)", "(m^3)", "(m)"))
                f.write("#\n")

        # Always append the new line of data
        time = self.timestep * self.delta_t
        volume = self.volume
        cg_x = self.com[0]
        cg_y = self.com[1]
        v_x = self.v_trans[0]
        v_y = self.v_trans[1]

        if self.rot_update == 'rot_mat':
            R = self.orientation[:2, :2]  # take 2D rotation part if orientation is 3x3
            theta_rad = np.arctan2(R[1, 0], R[0, 0])
        else:
            # Extract the w and z components from the orientation quaternion
            w = self.orientation[0]
            z = self.orientation[3]

            # The angle theta is 2 * arctan2(z, w)
            theta_rad = 2 * np.arctan2(z, w)

        # Convert to degrees for reporting
        theta_deg = np.degrees(theta_rad)

        force_x = self.force_int[0]
        force_y = self.force_int[1]
        moment_z = self.moment_int[2]
        h = self.h_min

        with open(file_name, "a") as f:
            f.write(f"{time:12.5e}  {cg_x:12.5e}  {cg_y:12.5e}  "
                    f"{v_x:12.5e}  {v_y:12.5e}  {theta_deg:12.5e}  "
                    f"{force_x:12.5e}  {force_y:12.5e}  {moment_z:12.5e}  {volume:12.5e}  {h:12.5e}\n")

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

    def calc_aitken(self, q_raw, q_prev, res_prev, relax_factor):
        """
        Generic Aitken Δ² relaxation factor update.
        Suitable for forces, moments, accelerations, or any vector residual.

        Args:
            q_raw       : Current unrelaxed vector (e.g. fluid force).
            q_prev      : Previous relaxed vector.
            res_prev    : Previous residual (q_raw_prev - q_prev_prev).
            relax_factor: Current relaxation factor.

        Returns:
            new_relax_factor (float): Updated and clamped Aitken relaxation factor.
        """

        # --- 1. STARTUP STABILIZATION ---
        if self.rb_iter <= 2:
            return relax_factor

        # --- 2. RAW AITKEN CALCULATION ---
        # Current residual
        res = q_raw - q_prev

        # Change in residual
        res_diff = res - res_prev

        # Compute denominator
        denom = np.dot(res_diff, res_diff)

        # Default to current factor if denominator is too small
        raw_relax = relax_factor

        # Update relaxation factor (only if safe)
        if denom > 1e-14:
            num = np.dot(res_prev, res_diff)
            raw_relax = - relax_factor * num / denom

        # --- 3. SMOOTHING (Slew Rate Limiter) ---
        # Limit the change in relaxation factor to avoid shocks (e.g. +/- 0.2 max change)
        max_change = 0.2

        change = raw_relax - relax_factor
        change = np.clip(change, -max_change, max_change)

        new_relax_factor = relax_factor + change

        # --- 4. SAFETY CLAMPING ---
        new_relax_factor = np.clip(new_relax_factor, 0.01, 0.99)

        return new_relax_factor
    
    def backup_fluent_log(self):
        file = join(self.dir_cfd, 'fluent.log')
        file_backup = join(self.dir_cfd, 'fluent_backup.log')
        if os.path.isfile(file_backup):
            os.remove(file_backup)
        if os.path.isfile(file):
            os.rename(file, file_backup)

    def print_rb_iteration_info(self, F_res, M_res):
        """
        Print residual info for rigid body iterations inside a coupling iteration.
        """

        info = f'   [Rigid body] {self.rb_iter:<18d}{F_res:<28.17e}{M_res:<28.17e}'
        tools.print_info(info, flush=True)

    def print_rb_header(self):
        """
        Print header for rigid body iterations inside a coupling iteration.
        """
        header = f'   {"RB Iteration":<18}{"Force residual":<28}{"Moment residual":<28}'
        tools.print_info(header, flush=True)

    def save_restart_rb_data(self):
        """Save rigid body motion state to a pickle file."""
        state = {
            'v_trans_prev': self.v_trans_prev,
            'v_trans': self.v_trans,
            'omega_prev': self.omega_prev,
            'omega': self.omega,
            'a_trans_prev': self.a_trans_prev,
            'a_trans': self.a_trans,
            'a_rot_prev': self.a_rot_prev,
            'a_rot': self.a_rot,
            'orientation_prev': self.orientation_prev,
            'orientation': self.orientation,
            'com': self.com,
            'com_prev': self.com_prev,
            'prev_patches': self.prev_patches,
            'new_patches': self.new_patches,
            'interface_rb': self.interface_rb,
            'force_pr_it': self.force_pr_it,
            'moment_pr_it': self.moment_pr_it
        }

        tmp = f'restart_rb_timestep{self.timestep}.pickle'
        file_name = join(self.dir_cfd, tmp)

        with open(file_name, 'wb') as f:
            pickle.dump(state, f)

    def load_restart_rb_data(self):
        """Load rigid body motion state from a pickle file."""
        if self.restart:
            tmp = f'restart_rb_timestep{self.timestep_start}.pickle'
        elif self.restart_rb_only != 0:
            tmp = f'restart_rb_timestep{self.restart_rb_only}.pickle'
        file_name = join(self.dir_cfd, tmp)

        if not os.path.exists(file_name):
            raise FileNotFoundError(f"Rigid body restart file not found: {file_name}")

        with open(file_name, 'rb') as f:
            state = pickle.load(f)

        self.v_trans_prev = state['v_trans_prev']
        self.v_trans = state['v_trans']
        self.omega_prev = state['omega_prev']
        self.omega = state['omega']
        self.a_trans_prev = state['a_trans_prev']
        self.a_trans = state['a_trans']
        self.a_rot_prev = state['a_rot_prev']
        self.a_rot = state['a_rot']
        self.com = state['com']
        self.com_prev = state['com_prev']
        self.prev_patches = state['prev_patches']
        self.new_patches = state['new_patches']
        self.force_pr_it = state['force_pr_it']
        self.moment_pr_it = state['moment_pr_it']
        if self.restart:
            self.interface_rb = state['interface_rb']
            self.orientation_prev = state['orientation_prev']
            self.orientation = state['orientation']

        tools.print_info('Rigid body restart data successfuly loaded.', layout='info')

# Helper functions
def quat_multiply(q1, q2):
    """Multiplies two quaternions."""
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    w = w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2
    x = w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2
    y = w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2
    z = w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2
    return np.array([w, x, y, z])

def quat_normalize(q):
    """Normalizes a quaternion to unit length to prevent numerical drift."""
    norm = np.linalg.norm(q)
    if norm == 0:
        return np.array([1.0, 0.0, 0.0, 0.0]) # Return identity quaternion
    return q / norm

def quat_from_angular_velocity(omega_vec, dt):
    """Creates a rotation quaternion from an angular velocity vector and timestep."""
    angle = np.linalg.norm(omega_vec) * dt
    if angle < 1e-12: # Avoid division by zero for very small rotations
        return np.array([1.0, 0.0, 0.0, 0.0])
    axis = omega_vec / (angle / dt)
    half_angle = angle / 2.0
    w = np.cos(half_angle)
    x, y, z = axis * np.sin(half_angle)
    return np.array([w, x, y, z])

def quat_rotate_vector_array(q, v_array):
    """Efficiently rotates an array of 3D vectors by a unit quaternion q."""
    q_w, q_vec = q[0], q[1:]
    # This is a fast, vectorized formula for quaternion rotation
    # Source: https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation?utm_source=chatgpt.com#Used_methods --> point 2
    return v_array + 2 * np.cross(q_vec, np.cross(q_vec, v_array) + q_w * v_array)

def skew(v):
    return np.array([[0, -v[2], v[1]],
                     [v[2], 0, -v[0]],
                     [-v[1], v[0], 0]])

def exp_map(x):
    if np.linalg.norm(x) < 1e-12:
        return np.identity(3)
    else:
        return np.identity(3) + (np.sin(np.linalg.norm(x)) / np.linalg.norm(x)) * skew(x) + (1 / 2) * (
                    (np.sin(np.linalg.norm(x)) ** 2) / ((np.linalg.norm(x) / 2) ** 2)) * (skew(x) @ skew(x))