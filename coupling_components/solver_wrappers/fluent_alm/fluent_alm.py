from coconut import data_structure
from coconut.coupling_components.solver_wrappers.solver_wrapper import SolverWrapper
from coconut import tools

import os
from os.path import join
import glob
import subprocess
import multiprocessing
import numpy as np
from scipy.interpolate import splprep, splev


def create(parameters):
    return SolverWrapperFluentALM(parameters)


class SolverWrapperFluentALM(SolverWrapper):
    # version specific parameters
    version = None  # Fluent product version, as from 2023R1 typically of the form 'xxxRx', set in subclass
    version_bis = None  # Fluent internal version, typically of the form 'x.x.0', set in subclass
    check_coupling_convergence_possible = True  # can solver check convergence after 1 iteration?

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
        self.cores = self.settings['cores']
        self.hosts_file = self.settings.get('hosts_file')
        self.case_file = self.settings['case_file']
        self.data_file = self.case_file.replace('.cas', '.dat', 1)
        if not os.path.exists(os.path.join(self.dir_cfd, self.case_file)):
            raise FileNotFoundError(f'Case file {self.case_file} not found in working directory {self.dir_cfd}')
        elif not os.path.exists(os.path.join(self.dir_cfd, self.data_file)):
            raise FileNotFoundError(f'Data file {self.data_file} not found in working directory {self.dir_cfd}')
        self.dimensions = 3  # hardcoded, UDF not yet adapted and tested for 2D-cases
        self.unsteady = self.settings['unsteady']
        self.flow_iterations = self.settings['flow_iterations']
        self.delta_t = self.settings['delta_t']
        self.timestep_start = self.settings['timestep_start']
        self.timestep = self.timestep_start
        self.save_results = self.settings.get('save_results', 1)
        self.save_restart = self.settings['save_restart']
        self.iteration = None
        self.fluent_process = None
        if len(self.settings['thread_names']) != 1:
            raise ValueError(f'Found {len(self.settings["thread_names"])} thread names, '
                             f'but solver supports only 1 object')
        self.thread_ids = {}  # thread IDs corresponding to thread names
        # Since threads not physically present, give them custom id's
        for i, thread_name in enumerate(self.settings['thread_names']):
            self.thread_ids[thread_name] = i
        self.model_part_thread_ids = {}  # thread IDs corresponding to ModelParts
        self.model = None

        self.alm_settings = self.settings['ALM']
        self.nyarnpoints = None
        self.yarn_diameter = self.alm_settings['yarn_diameter']
        self.g_eps = self.alm_settings['g_eps']  # shape parameter of force smearing kernel (2D Gaussian)
        self.n_circ_s = self.alm_settings.get('n_circ_s', 5)  # # circular sampling points for axial velocity sampling

    @tools.time_initialize
    def initialize(self):
        super().initialize()

        # prepare Fluent journal
        journal = f'alm.jou'
        unsteady = '#t' if self.unsteady else '#f'
        check_coupling_convergence = '#t' if self.check_coupling_convergence else '#f'
        with open(join(self.dir_src, journal)) as infile:
            with open(join(self.dir_cfd, journal), 'w') as outfile:
                for line in infile:
                    line = line.replace('|CASE|', join(self.dir_cfd, self.case_file))
                    line = line.replace('|UNSTEADY|', unsteady)
                    line = line.replace('|FLOW_ITERATIONS|', str(self.flow_iterations))
                    line = line.replace('|CHECK_COUPLING_CONVERGENCE|', check_coupling_convergence)
                    line = line.replace('|DELTA_T|', str(self.delta_t))
                    line = line.replace('|TIMESTEP_START|', str(self.timestep_start))
                    line = line.replace('|END_OF_TIMESTEP_COMMANDS|', self.settings.get('end_of_timestep_commands',
                                                                                        '\n'))
                    outfile.write(line)

        # read in coordinate file to determine number of points in the yarn
        # format of file: nyarnpoints+1 rows and self.dimensions columns with x-, y- (and z-) coordinates in m;
        # first row: number of yarn points (nyarnpoints)
        file_name = join(self.dir_cfd, f'coordinates_timestep0.dat')
        with open(file_name, 'r') as f:
            self.nyarnpoints = int(f.readline().strip())
        data = np.loadtxt(file_name, skiprows=1)
        if data.shape[1] != self.dimensions:
            raise ValueError('Given dimension does not match coordinates')
        if data.shape[0] != self.nyarnpoints:
            raise ValueError(f'Declared number of points ({self.nyarnpoints}) does not match number of rows '
                             f'({data.shape[0]})')

        # prepare Fluent UDF
        udf = f'alm.c'
        with open(join(self.dir_src, udf)) as infile:
            with open(join(self.dir_cfd, udf), 'w') as outfile:
                for line in infile:
                    line = line.replace('|NYARNPOINTS|', str(self.nyarnpoints))
                    line = line.replace('|YARN_DIAMETER|', str(self.yarn_diameter))
                    line = line.replace('|G_EPS|', str(self.g_eps))
                    line = line.replace('|N_CIRC_S|', str(self.n_circ_s))
                    line = line.replace('|DELTA_T|', str(self.delta_t))
                    line = line.replace('|UNSTEADY|', str(int(self.unsteady)))
                    outfile.write(line)

        # check number of cores
        if self.hosts_file is not None:
            with open(join(self.dir_cfd, self.hosts_file)) as fp:
                max_cores = len(fp.readlines())
        else:
            max_cores = multiprocessing.cpu_count()
        if self.cores < 1 or self.cores > max_cores:
            tools.print_info(f'Number of cores incorrect, changed from {self.cores} to {max_cores}', layout='warning')
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
        self.fluent_process = subprocess.Popen(cmd, executable='/bin/bash', shell=True, cwd=self.dir_cfd, env=self.env)

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
                            raise ValueError('Dimension in JSON does not match Fluent')
                    check = 2
                elif check == 2 and 'Time' in line:
                    if 'Steady' in line and self.unsteady:
                        raise ValueError('Unsteady in JSON does not match steady Fluent')
                    elif 'Unsteady' in line and not self.unsteady:
                        raise ValueError('Steady in JSON does not match unsteady Fluent')

        if os.path.isfile(join(self.dir_cfd, 'log')):
            os.unlink(join(self.dir_cfd, 'log'))  # delete log file (fluent.log is sufficient)

        # remove "report.sum" because the batch options to overwrite report files and case files conflict in some
        # versions of Fluent (2023R1)
        os.unlink(report)

        # create Model
        self.model = data_structure.Model()

        # create input ModelParts (nodes)
        for item in (self.settings['interface_input']):
            mp_name = item['model_part']

            # get face thread ID that corresponds to ModelPart
            for thread_name in self.thread_ids:
                if thread_name in mp_name:
                    self.model_part_thread_ids[mp_name] = self.thread_ids[thread_name]
            if mp_name not in self.model_part_thread_ids:
                raise AttributeError('Could not find thread name corresponding ' +
                                     f'to ModelPart {mp_name}')

            # get node coordinates and ids (data file is still in memory), assume that nodes are sorted by connectivity
            coords_tmp = np.zeros((data.shape[0], 3)) * 0.
            coords_tmp[:, :self.dimensions] = data  # add column z if 2D
            ids_tmp = np.arange(self.nyarnpoints).astype(int)  # create node id's

            # sort and remove doubles (in fact obsolete)
            args = np.unique(ids_tmp, return_index=True)[1].tolist()
            x0 = coords_tmp[args, 0]
            y0 = coords_tmp[args, 1]
            z0 = coords_tmp[args, 2]
            ids = ids_tmp[args]

            # create ModelPart
            self.model.create_model_part(mp_name, x0, y0, z0, ids)

        # create output ModelParts (faces)
        for item in (self.settings['interface_output']):
            mp_name = item['model_part']

            # get face thread ID that corresponds to ModelPart
            for thread_name in self.thread_ids:
                if thread_name in mp_name:
                    self.model_part_thread_ids[mp_name] = self.thread_ids[thread_name]
            if mp_name not in self.model_part_thread_ids:
                raise AttributeError('Could not find thread name corresponding ' +
                                     f'to ModelPart {mp_name}')

            # get node coordinates and ids (data file is still in memory), assume that nodes are sorted by connectivity
            coords_tmp = np.zeros((data.shape[0], 3)) * 0.
            coords_tmp[:, :self.dimensions] = data  # add column z if 2D
            ids_tmp = np.arange(self.nyarnpoints).astype(int)  # create node id's

            # sort and remove doubles (in fact obsolete)
            args = np.unique(ids_tmp, return_index=True)[1].tolist()
            x0 = coords_tmp[args, 0]
            y0 = coords_tmp[args, 1]
            z0 = coords_tmp[args, 2]
            ids = ids_tmp[args]

            # create ModelPart
            self.model.create_model_part(mp_name, x0, y0, z0, ids)

        # create Interfaces
        self.interface_input = data_structure.Interface(self.settings['interface_input'], self.model)
        self.interface_output = data_structure.Interface(self.settings['interface_output'], self.model)

        # write node positions and let Fluent know that ModelParts are created
        if self.timestep_start == 0:
            self.write_node_positions()
        else:  # restart: mandatory to keep last coordinates_update-file
            if not os.path.exists(join(self.dir_cfd, f'coordinates_update_timestep{self.timestep_start}.dat')):
                raise FileNotFoundError(f'Could not find file coordinates_update_timestep{self.timestep_start}.dat')
        self.coco_messages.send_message('model_parts_created')

    def initialize_solution_step(self):
        super().initialize_solution_step()

        self.iteration = 0
        self.timestep += 1

        self.coco_messages.send_message('next')
        self.coco_messages.wait_message('next_ready')

    @tools.time_solve_solution_step
    def solve_solution_step(self, interface_input):
        self.iteration += 1

        # store incoming displacements
        self.interface_input.set_interface_data(interface_input.get_interface_data())

        # write interface data
        self.write_node_positions()

        # copy input data for debugging
        if self.debug:
            for dct in self.interface_input.parameters:
                mp_name = dct['model_part']
                thread_id = self.model_part_thread_ids[mp_name]
                src = f'coordinates_update_timestep{self.timestep}.dat'
                dst = f'coordinates_update_timestep{self.timestep}_Iter{self.iteration}.dat'
                cmd = f'cp {join(self.dir_cfd, src)} {join(self.dir_cfd, dst)}'
                os.system(cmd)

        # let Fluent run, wait for data
        self.coco_messages.send_message('continue')
        self.coco_messages.wait_message('continue_ready')

        if self.check_coupling_convergence:
            # check if Fluent converged after 1 iteration
            self.coupling_convergence = self.coco_messages.check_message('solver_converged')
            if self.print_coupling_convergence and self.coupling_convergence:
                tools.print_info(f'{self.__class__.__name__} converged')

        # read data from Fluent
        for dct in self.interface_output.parameters:
            mp_name = dct['model_part']

            # read in datafile
            tmp = f'traction_timestep{self.timestep}.dat'
            file_name = join(self.dir_cfd, tmp)
            data = np.loadtxt(file_name, skiprows=1)
            if data.shape[1] != self.dimensions + 1:
                raise ValueError('Given dimension does not match coordinates')

            # copy output data for debugging
            if self.debug:
                dst = f'traction_timestep{self.timestep}_it{self.iteration}.dat'
                cmd = f'cp {file_name} {join(self.dir_cfd, dst)}'
                os.system(cmd)

            # get face coordinates and ids
            traction_tmp = np.zeros((data.shape[0], 3)) * 0.
            traction_tmp[:, :self.dimensions] = data[:, :-1]
            ids_tmp = data[:, -1]

            # sort and remove doubles
            args = np.unique(ids_tmp, return_index=True)[1].tolist()
            traction = traction_tmp[args, :]
            ids = ids_tmp[args]

            # store pressure and traction in Nodes
            model_part = self.model.get_model_part(mp_name)
            if ids.size != model_part.size:
                raise ValueError('Size of data does not match size of ModelPart')
            if not np.all(ids == model_part.id):
                raise ValueError('IDs of data do not match ModelPart IDs')

            self.interface_output.set_variable_data(mp_name, 'traction', traction)

        # return interface_output object
        return self.interface_output

    def finalize_solution_step(self):
        super().finalize_solution_step()

    @tools.time_save
    def output_solution_step(self):
        super().output_solution_step()

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
                try:
                    os.remove(join(self.dir_cfd, f'coordinates_update_timestep{self.timestep + self.save_restart}.dat'))
                except OSError:
                    pass

    def finalize(self):
        super().finalize()
        self.coco_messages.send_message('stop')
        self.coco_messages.wait_message('stop_ready')
        self.fluent_process.wait()

        # remove unnecessary files
        self.remove_dat_files(self.timestep)

        # delete .trn files
        for path in glob.glob(join(self.dir_cfd, '*.trn')):
            os.remove(path)

    def remove_dat_files(self, timestep):
        if not self.debug:
            try:
                os.remove(join(self.dir_cfd, f'traction_timestep{timestep}.dat'))
            except OSError:
                pass

    def check_software(self):
        # Fluent version: see set_fluent_version
        result = subprocess.run(['fluent', '-r'], stdout=subprocess.PIPE, env=self.env)
        if self.version_bis not in str(result.stdout):
            raise RuntimeError(f'ANSYS Fluent version {self.version} ({self.version_bis}) is required. Check if '
                               f'the solver load commands for the "machine_name" are correct in solver_modules.py.')

    def write_node_positions(self):
        for dct in self.interface_input.parameters:
            mp_name = dct['model_part']
            model_part = self.model.get_model_part(mp_name)
            if self.timestep == 0:
                x = model_part.x0
                y = model_part.y0
                z = model_part.z0
            else:
                displacement = self.interface_input.get_variable_data(mp_name, 'displacement')
                x = model_part.x0 + displacement[:, 0]
                y = model_part.y0 + displacement[:, 1]
                z = model_part.z0 + displacement[:, 2]
            tck_u = splprep([x, y, z], s=0, k=3)
            t = np.array(splev(tck_u[1], tck_u[0], der=1)).T
            if self.dimensions == 2:
                data = np.rec.fromarrays([x, y, t[:, 0], t[:, 1]])
                fmt = '%27.17e%27.17e%27.17e%27.17e'
            else:
                data = np.rec.fromarrays([x, y, z, t[:, 0], t[:, 1], t[:, 2]])
                fmt = '%27.17e%27.17e%27.17e%27.17e%27.17e%27.17e'
            tmp = f'coordinates_update_timestep{self.timestep}.dat'
            file_name = join(self.dir_cfd, tmp)
            np.savetxt(file_name, data, fmt=fmt, header=f'{model_part.size}', comments='')

    def get_coordinates(self):
        """  # TODO: rewrite this
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

            # read in datafile
            tmp = f'nodes_timestep{self.timestep}.dat'
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

        # no ids and coordinates for output ModelParts as this has by definition the same coordinates as the input MPs

        return coord_data

    def backup_fluent_log(self):
        file = join(self.dir_cfd, 'fluent.log')
        file_backup = join(self.dir_cfd, 'fluent_backup.log')
        if os.path.isfile(file_backup):
            os.remove(file_backup)
        if os.path.isfile(file):
            os.rename(file, file_backup)
