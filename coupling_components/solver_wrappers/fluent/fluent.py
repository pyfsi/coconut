from coconut import data_structure
from coconut.coupling_components.component import Component
from coconut import tools

import os
from os.path import join
import glob
import subprocess
import multiprocessing
import time
import numpy as np
import sys
import hashlib
from getpass import getuser
import shutil

# TODO: issue: id and hash shadow built-in names


def create(parameters):
    return SolverWrapperFluent(parameters)


class SolverWrapperFluent(Component):
    # version specific parameters
    version = None  # Fluent product version, as from 2019R1 typically of the form 'xxxRx', set in sub-class
    version_bis = None  # Fluent internal version, typically of the form 'x.x.0', set in sub-class

    @tools.time_initialize
    def __init__(self, parameters):
        super().__init__()

        if self.version is None or self.version_bis is None:
            raise NotImplementedError(
                'Base class method called, class variable version and version_bis need to be set in the derived class')

        # set parameters
        self.settings = parameters['settings']
        self.dir_cfd = join(os.getcwd(), self.settings['working_directory'])
        self.env = None  # environment in which correct version of Fluent software is available, set in sub-class
        self.remove_all_messages()
        self.backup_fluent_log()
        self.dir_src = os.path.realpath(os.path.dirname(__file__))
        self.tmp_directory_name = f'coconut_{getuser()}_{os.getpid()}_fluent'  # dir in /tmp for host-node communication
        self.cores = self.settings['cores']
        self.hostfile = self.settings.get('hostfile', False)
        if self.hostfile:
            with open(self.hostfile) as fp:
                max_cores = len(fp.readlines())
        else:
            max_cores = multiprocessing.cpu_count()
        if self.cores < 1 or self.cores > max_cores:
            self.cores = max_cores
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
        self.model = None
        self.interface_input = None
        self.interface_output = None

        # time
        self.init_time = self.init_time
        self.run_time = 0.0

        # debug
        self.debug = False  # set on True to save copy of input and output files in every iteration

    @tools.time_initialize
    def initialize(self):
        super().initialize()

        # prepare Fluent journal
        journal = f'v{self.version}.jou'
        thread_names_str = ''
        for thread_name in self.thread_ids:
            thread_names_str += ' "' + thread_name + '"'
        unsteady = '#f'
        if self.unsteady:
            unsteady = '#t'
        multiphase = '#f'
        if self.multiphase:
            multiphase = '#t'
        with open(join(self.dir_src, journal)) as infile:
            with open(join(self.dir_cfd, journal), 'w') as outfile:
                for line in infile:
                    line = line.replace('|CASE|', join(self.dir_cfd, self.case_file))
                    line = line.replace('|THREAD_NAMES|', thread_names_str)
                    line = line.replace('|UNSTEADY|', unsteady)
                    line = line.replace('|MULTIPHASE|', multiphase)
                    line = line.replace('|FLOW_ITERATIONS|', str(self.flow_iterations))
                    line = line.replace('|DELTA_T|', str(self.delta_t))
                    line = line.replace('|TIMESTEP_START|', str(self.timestep_start))
                    line = line.replace('|END_OF_TIMESTEP_COMMANDS|', self.settings.get('end_of_timestep_commands',
                                                                                        '\n'))
                    outfile.write(line)

        # prepare Fluent UDF
        udf = f'v{self.version}.c'
        shutil.rmtree(join('/tmp', self.tmp_directory_name), ignore_errors=True)
        os.mkdir(join('/tmp', self.tmp_directory_name))
        with open(join(self.dir_src, udf)) as infile:
            with open(join(self.dir_cfd, udf), 'w') as outfile:
                for line in infile:
                    line = line.replace('|MAX_NODES_PER_FACE|', str(self.mnpf))
                    line = line.replace('|TMP_DIRECTORY_NAME|', self.tmp_directory_name)
                    outfile.write(line)

        # start Fluent with journal
        log = join(self.dir_cfd, 'fluent.log')
        cmd1 = f'fluent -r{self.version_bis} {self.dimensions}ddp '
        cmd2 = f'-t{self.cores} -i {journal}'

        if self.hostfile:
            cmd1 += f' -cnf={self.hostfile} -ssh '
        if self.settings['fluent_gui']:
            cmd = cmd1 + cmd2
        else:
            cmd = cmd1 + '-gu ' + cmd2 + f' >> {log} 2>&1'  # TODO: does log work well?
            # cmd = cmd1 + '-gu ' + cmd2 + f' 2> >(tee -a {log}) 1>> {log}'  # TODO: specify what this option does?
        self.fluent_process = subprocess.Popen(cmd, executable='/bin/bash',
                                               shell=True, cwd=self.dir_cfd, env=self.env)

        # get general simulation info from  fluent.log and report.sum
        self.wait_message('case_info_exported')

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
                    elif'Numerics' in line:
                        if self.multiphase:
                            raise ValueError('Multiphase in JSON does not match singlephase Fluent')
                        break

        if os.path.isfile(join(self.dir_cfd, 'log')):
            os.unlink(join(self.dir_cfd, 'log'))  # delete log file (fluent.log is sufficient)

        # get surface thread ID's from report.sum and write them to bcs.txt
        check = 0
        info = []
        with open(report, 'r') as file:
            for line in file:
                if check == 3 and line.islower():
                    name, id, _ = line.strip().split()
                    if name in self.thread_ids:
                        info.append(' '.join((name, id)))
                        self.thread_ids[name] = id

                if check == 3 and not line.islower():
                    break
                if check == 2:  # skip 1 line
                    check = 3
                if 'name' in line and check == 1:
                    check = 2
                if 'Boundary Conditions' in line:
                    check = 1
        with open(join(self.dir_cfd, 'bcs.txt'), 'w') as file:
            file.write(str(len(info)) + '\n')
            for line in info:
                file.write(line + '\n')
        self.send_message('thread_ids_written_to_file')

        # import node and face information if no restart
        if self.timestep_start == 0:
            self.wait_message('nodes_and_faces_stored')

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

            # read in datafile
            thread_id = self.model_part_thread_ids[mp_name]
            file_name = join(self.dir_cfd, f'nodes_timestep0_thread{thread_id}.dat')
            data = np.loadtxt(file_name, skiprows=1)
            if data.shape[1] != self.dimensions + 1:
                raise ValueError('Given dimension does not match coordinates')

            # get node coordinates and ids
            coords_tmp = np.zeros((data.shape[0], 3)) * 0.
            coords_tmp[:, :self.dimensions] = data[:, :-1]  # add column z if 2D
            ids_tmp = data[:, -1].astype(int)  # array is flattened

            # sort and remove doubles
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
            # read in datafile
            thread_id = self.model_part_thread_ids[mp_name]
            file_name = join(self.dir_cfd, f'faces_timestep0_thread{thread_id}.dat')
            data = np.loadtxt(file_name, skiprows=1)
            if data.shape[1] != self.dimensions + self.mnpf:
                raise ValueError(f'Given dimension does not match coordinates')

            # get face coordinates and ids
            coords_tmp = np.zeros((data.shape[0], 3)) * 0.
            coords_tmp[:, :self.dimensions] = data[:, :-self.mnpf]  # add column z if 2D
            ids_tmp = self.get_unique_face_ids(data[:, -self.mnpf:])

            # sort and remove doubles
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

    def initialize_solution_step(self):
        super().initialize_solution_step()

        self.iteration = 0
        self.timestep += 1

        self.send_message('next')
        self.wait_message('next_ready')

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
                src = f'nodes_update_timestep{self.timestep}_thread{thread_id}.dat'
                dst = f'nodes_update_timestep{self.timestep}_thread{thread_id}_Iter{self.iteration}.dat'
                cmd = f'cp {join(self.dir_cfd, src)} {join(self.dir_cfd, dst)}'
                os.system(cmd)

        # let Fluent run, wait for data
        self.send_message('continue')
        self.wait_message('continue_ready')

        # read data from Fluent
        for dct in self.interface_output.parameters:
            mp_name = dct['model_part']

            # read in datafile
            thread_id = self.model_part_thread_ids[mp_name]
            tmp = f'pressure_traction_timestep{self.timestep}_thread{thread_id}.dat'
            file_name = join(self.dir_cfd, tmp)
            data = np.loadtxt(file_name, skiprows=1)
            if data.shape[1] != self.dimensions + 1 + self.mnpf:
                raise ValueError('Given dimension does not match coordinates')

            # copy output data for debugging
            if self.debug:  # TODO: Iter --> iter everywhere?
                dst = f'pressure_traction_timestep{self.timestep}_thread{thread_id}_Iter{self.iteration}.dat'
                cmd = f'cp {file_name} {join(self.dir_cfd, dst)}'
                os.system(cmd)

            # get face coordinates and ids
            traction_tmp = np.zeros((data.shape[0], 3)) * 0.
            traction_tmp[:, :self.dimensions] = data[:, :-1 - self.mnpf]
            pressure_tmp = data[:, self.dimensions].reshape(-1, 1)
            ids_tmp = self.get_unique_face_ids(data[:, -self.mnpf:])

            # sort and remove doubles
            args = np.unique(ids_tmp, return_index=True)[1].tolist()
            traction = traction_tmp[args, :]
            pressure = pressure_tmp[args]
            ids = ids_tmp[args]

            # store pressure and traction in Nodes
            model_part = self.model.get_model_part(mp_name)
            if ids.size != model_part.size:
                raise ValueError('Size of data does not match size of ModelPart')
            if not np.all(ids == model_part.id):
                raise ValueError('IDs of data do not match ModelPart IDs')

            self.interface_output.set_variable_data(mp_name, 'traction', traction)
            self.interface_output.set_variable_data(mp_name, 'pressure', pressure)

        # return interface_output object
        return self.interface_output

    def finalize_solution_step(self):
        super().finalize_solution_step()

        # save if required
        if (self.save_results != 0 and self.timestep % self.save_results == 0) \
                or (self.save_restart != 0 and self.timestep % self.save_restart == 0):
            self.send_message('save')
            self.wait_message('save_ready')

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
        shutil.rmtree(join('/tmp', self.tmp_directory_name), ignore_errors=True)
        self.send_message('stop')
        self.wait_message('stop_ready')
        self.fluent_process.wait()

        # remove unnecessary files
        self.remove_dat_files(self.timestep)

        # delete .trn files
        for path in glob.glob(join(self.dir_cfd, '*.trn')):
            os.remove(path)

    def remove_dat_files(self, timestep):
        if not self.debug:
            for thread_id in self.thread_ids.values():
                try:
                    os.remove(join(self.dir_cfd, f'nodes_update_timestep{timestep}_thread{thread_id}.dat'))
                    os.remove(join(self.dir_cfd, f'pressure_traction_timestep{timestep}_thread{thread_id}.dat'))
                except OSError:
                    pass

    def get_interface_input(self):
        return self.interface_input

    def get_interface_output(self):
        return self.interface_output

    def check_software(self):
        # Python version: 3.6 or higher
        if sys.version_info < (3, 6):
            raise RuntimeError('Python version 3.6 or higher required.')

        # Fluent version: see set_fluent_version
        result = subprocess.run(['fluent', '-r'], stdout=subprocess.PIPE, env=self.env)
        if self.version_bis not in str(result.stdout):
            raise RuntimeError(f'ANSYS Fluent version {self.version} ({self.version_bis}) is required')

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
        """
        data = data.astype(int)
        # ids = np.zeros(data.shape[0], dtype='U256')  # array is flattened
        ids = np.zeros(data.shape[0], dtype=int)
        for i in range(ids.size):
            tmp = np.unique(data[i, :])
            if tmp[0] == -1:
                tmp = tmp[1:]
            unique_string = '-'.join(tuple(tmp.astype(str)))
            hash = hashlib.sha1(str.encode(unique_string))
            ids[i] = int(hash.hexdigest(), 16) % (10 ** 16)
        return ids

    def write_node_positions(self):
        for dct in self.interface_input.parameters:
            mp_name = dct['model_part']
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
        self.send_message('store_grid')
        self.wait_message('store_grid_ready')

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
            data = np.loadtxt(join(self.dir_cfd, tmp), skiprows=1)
            if data.shape[1] != self.dimensions + self.mnpf:
                raise ValueError(f'Given dimension does not match coordinates')

            # get face coordinates and ids
            coords_tmp = np.zeros((data.shape[0], 3)) * 0.
            coords_tmp[:, :self.dimensions] = data[:, :-self.mnpf]  # add column z if 2D
            ids_tmp = self.get_unique_face_ids(data[:, -self.mnpf:])

            # sort and remove doubles
            args = np.unique(ids_tmp, return_index=True)[1].tolist()
            coord_data[mp_name]['ids'] = ids_tmp[args]
            coord_data[mp_name]['coords'] = coords_tmp[args, :]

        return coord_data

    def send_message(self, message):
        file = join(self.dir_cfd, message + '.coco')
        open(file, 'w').close()
        return

    def wait_message(self, message):
        file = join(self.dir_cfd, message + '.coco')
        while not os.path.isfile(file):
            time.sleep(0.01)
        os.remove(file)
        return

    def check_message(self, message):
        file = join(self.dir_cfd, message + '.coco')
        if os.path.isfile(file):
            os.remove(file)
            return True
        return False

    def remove_all_messages(self):
        for file_name in os.listdir(self.dir_cfd):
            if file_name.endswith('.coco'):
                file = join(self.dir_cfd, file_name)
                os.remove(file)

    def backup_fluent_log(self):
        file = join(self.dir_cfd, 'fluent.log')
        file_backup = join(self.dir_cfd, 'fluent_backup.log')
        if os.path.isfile(file_backup):
            os.remove(file_backup)
        if os.path.isfile(file):
            os.rename(file, file_backup)
