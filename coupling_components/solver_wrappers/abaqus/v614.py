from coconut import data_structure
from coconut.coupling_components.component import Component
from coconut import tools

import os
from os.path import join
import subprocess
import shutil
import sys
import time
import numpy as np
import re
import textwrap
from pathlib import Path


def create(parameters):
    return SolverWrapperAbaqus614(parameters)


class SolverWrapperAbaqus614(Component):
    version = '6.14'

    def __init__(self, parameters):
        super().__init__()

        # set parameters
        self.check_software()
        self.settings = parameters['settings']
        self.dir_csm = join(os.getcwd(), self.settings['working_directory'])  # *** alternative for getcwd?
        path_src = os.path.realpath(os.path.dirname(__file__))
        self.logfile = 'abaqus.log'

        self.cores = self.settings['cores']  # number of CPUs Abaqus has to use
        self.dimensions = self.settings['dimensions']
        self.array_size = self.settings['arraysize']
        self.delta_t = self.settings['delta_t']
        self.timestep_start = self.settings['timestep_start']
        self.surfaceIDs = self.settings['surfaceIDs']
        self.n_surfaces = len(self.surfaceIDs)
        self.mp_mode = self.settings['mp_mode']
        self.input_file = self.settings['input_file']
        self.timestep = self.timestep_start
        self.iteration = None
        self.model_part_surface_ids = {}  # surface IDs corresponding to ModelParts

        self.subcycling = self.settings.get('subcycling', False)  # value from parameters or False if not present
        if self.subcycling:
            self.min_inc = self.settings['min_inc']
            self.initial_inc = self.settings['initial_inc']
            self.max_num_inc = self.settings['max_num_inc']
            self.max_inc = self.settings['max_inc']
            self.ramp = (0, 1)[self.settings['ramp']]  # 0 or 1 required to substitute in user-subroutines (FORTRAN)
        else:
            self.ramp = 0
            self.max_num_inc = 1

        # prepare abaqus_v6.env
        hostname_replace = ''
        if self.mp_mode == 'MPI' and 'AbaqusHosts.txt' in os.listdir(self.dir_csm):
            with open(join(self.dir_csm, 'AbaqusHosts.txt'), 'r') as host_file:
                host_names = host_file.read().split()
            hostname_replace = str([[hostname, host_names.count(hostname)] for hostname in set(host_names)])
        with open(join(path_src, 'abaqus_v6.env'), 'r') as infile:
            with open(join(self.dir_csm, 'abaqus_v6.env'), 'w') as outfile:
                for line in infile:
                    line = line.replace('|HOSTNAME|', hostname_replace) if self.mp_mode == 'MPI' \
                        else (line, '')['|HOSTNAME|' in line]  # replace |HOSTNAME| if MPI else remove line
                    line = line.replace('|MP_MODE|', self.mp_mode)
                    line = line.replace('|PID|', str(os.getpid()))
                    if '|' in line:
                        raise ValueError(f'The following line in abaqus_v6.env still contains a \'|\' after '
                                         f'substitution: \n \t{line} \n Probably a parameter was not substituted')
                    outfile.write(line)

        # create start and restart file (each time step > 1 is a restart of the previous converged time step)
        self.write_start_and_restart_inp(join(self.dir_csm, self.input_file), self.dir_csm + '/CSM_Time0.inp',
                                         self.dir_csm + '/CSM_Restart.inp')

        # prepare Abaqus USRInit.f
        usr = 'USRInit.f'
        with open(join(path_src, usr), 'r') as infile:
            with open(join(self.dir_csm, 'usrInit.f'), 'w') as outfile:
                for line in infile:
                    line = line.replace('|dimension|', str(self.dimensions))
                    line = line.replace('|surfaces|', str(self.n_surfaces))
                    line = line.replace('|cpus|', str(self.cores))

                    # if PWD is too long then FORTRAN code can not compile so this needs special treatment
                    line = self.replace_fortran(line, '|PWD|', os.path.abspath(os.getcwd()))
                    line = self.replace_fortran(line, '|CSM_dir|', self.settings['working_directory'])
                    if '|' in line:
                        raise ValueError(f'The following line in USRInit.f still contains a \'|\' after substitution: '
                                         f'\n \t{line} \nProbably a parameter was not substituted')
                    outfile.write(line)

        # compile Abaqus USRInit.f in library libusr
        path_libusr = join(self.dir_csm, 'libusr/')
        os.system('rm -rf ' + path_libusr)
        os.system('mkdir ' + path_libusr)
        cmd = f'abaqus make library=usrInit.f directory={path_libusr} >> {self.logfile} 2>&1'
        self.print_log(f'### Compilation of usrInit.f ###')
        subprocess.run(cmd, shell=True, cwd=self.dir_csm, executable='/bin/bash')

        # get load points from usrInit.f at timestep_start
        self.print_log(f'\n### Get load integration points using usrInit.f ###')
        if self.timestep_start == 0:
            cmd1 = f'rm -f CSM_Time{self.timestep_start}Surface*Faces.dat ' \
                f'CSM_Time{self.timestep_start}Surface*FacesBis.dat'
            # the output files will have a name with a higher time step  ('job=') than the input file ('input=')
            cmd2 = f'abaqus job=CSM_Time{self.timestep_start + 1} input=CSM_Time{self.timestep_start} ' \
                f'cpus=1 output_precision=full interactive >> {self.logfile} 2>&1'
            commands = cmd1 + '; ' + cmd2
            subprocess.run(commands, shell=True, cwd=self.dir_csm, executable='/bin/bash')
        else:  # restart: this only used for checks
            cmd1 = f'rm -f CSM_Time{self.timestep_start}Surface*Faces.dat ' \
                f'CSM_Time{self.timestep_start}Surface*FacesBis.dat'
            cmd2 = f'abaqus job=CSM_Time{self.timestep_start + 1} oldjob=CSM_Time{self.timestep_start} ' \
                f'input=CSM_Restart cpus=1 output_precision=full interactive >> {self.logfile} 2>&1'
            commands = cmd1 + '; ' + cmd2
            subprocess.run(commands, shell=True, cwd=self.dir_csm, executable='/bin/bash')

        # prepare GetOutput.cpp
        get_output = 'GetOutput.cpp'
        temp_str = ''
        for j in range(0, self.n_surfaces - 1):
            temp_str += f'\"{self.surfaceIDs[j]}\", '
        temp_str += f'\"{self.surfaceIDs[self.n_surfaces-1]}\"'

        with open(join(path_src, get_output), 'r') as infile:
            with open(join(self.dir_csm, get_output), 'w') as outfile:
                for line in infile:
                    line = line.replace('|surfaces|', str(self.n_surfaces))
                    line = line.replace('|surfaceIDs|', temp_str)
                    line = line.replace('|dimension|', str(self.dimensions))
                    if '|' in line:
                        raise ValueError(f'The following line in GetOutput.cpp still contains a \'|\' after '
                                         f'substitution: \n \t{line} \n Probably a parameter was not substituted')
                    outfile.write(line)

        # compile GetOutput.cpp
        self.print_log(f'\n### Compilation of GetOutput.cpp ###')
        cmd = f'abaqus make job=GetOutput user=GetOutput.cpp >> {self.logfile} 2>&1'
        subprocess.run(cmd, shell=True, cwd=self.dir_csm, executable='/bin/bash')

        # get node positions (not load points) at timestep_start (0 is an argument to GetOutput.exe)
        self.print_log(f'\n### Get geometrical node positions using GetOutput ###')
        cmd = f'abaqus ./GetOutput.exe CSM_Time{self.timestep_start + 1} 0 >> {self.logfile} 2>&1'
        subprocess.run(cmd, shell=True, cwd=self.dir_csm, executable='/bin/bash')

        for i in range(0, self.n_surfaces):
            path_output = join(self.dir_csm, f'CSM_Time{self.timestep_start + 1}Surface{i}Output.dat')
            path_nodes = join(self.dir_csm, f'CSM_Time{self.timestep_start}Surface{i}Nodes.dat')
            shutil.move(path_output, path_nodes)

            # create elements file per surface
            face_file = os.path.join(self.dir_csm, f'CSM_Time{self.timestep_start}Surface{i}Cpu0Faces.dat')
            output_file = os.path.join(self.dir_csm, f'CSM_Time{self.timestep_start}Surface{i}Elements.dat')
            self.make_elements(face_file, output_file)

        # prepare Abaqus USR.f
        usr = 'USR.f'
        with open(join(path_src, usr), 'r') as infile:
            with open(join(self.dir_csm, 'usr.f'), 'w') as outfile:
                for line in infile:
                    line = line.replace('|dimension|', str(self.dimensions))
                    line = line.replace('|arraySize|', str(self.array_size))
                    line = line.replace('|surfaces|', str(self.n_surfaces))
                    line = line.replace('|cpus|', str(self.cores))
                    line = line.replace('|ramp|', str(self.ramp))
                    line = line.replace('|deltaT|', str(self.delta_t))

                    # if PWD is too long then FORTRAN code cannot compile so this needs special treatment
                    line = self.replace_fortran(line, '|PWD|', os.path.abspath(os.getcwd()))
                    line = self.replace_fortran(line, '|CSM_dir|', self.settings['working_directory'])
                    if '|' in line:
                        raise ValueError(f'The following line in USR.f still contains a \'|\' after substitution: '
                                         f'\n \t{line} \n Probably a parameter was not substituted')
                    outfile.write(line)

        # compile Abaqus USR.f
        self.print_log(f'\n### Compilation of usr.f ###')
        os.system('rm -r ' + path_libusr)  # remove libusr containing compiled USRInit.f
        os.system('mkdir ' + path_libusr)
        cmd = f'abaqus make library=usr.f directory={path_libusr} >> {self.logfile} 2>&1'
        subprocess.run(cmd, shell=True, cwd=self.dir_csm, executable='/bin/bash')

        # create Model
        self.model = data_structure.Model()

        # create input ModelParts (load points)
        for item in (self.settings['interface_input']):
            mp_name = item['model_part']

            for i, surfaceID in enumerate(self.surfaceIDs):  # identify surfaceID corresponding to ModelPart
                if surfaceID in mp_name:
                    self.model_part_surface_ids[mp_name] = i
                    break
            if mp_name not in self.model_part_surface_ids:
                raise AttributeError(f'Could not identify surfaceID corresponding to ModelPart {mp_name}')
            mp_id = self.model_part_surface_ids[mp_name]

            # read in elements file
            elem0_file = join(self.dir_csm, f'CSM_Time0Surface{mp_id}Elements.dat')
            elements0 = np.loadtxt(elem0_file)
            n_elem = int(elements0[0, 0])  # elements first item on line 1 contains number of elements
            n_lp = int(elements0[0, 1])  # elements second item on line 1 contains number of total load points
            if elements0.shape[0] - 1 != int(n_elem):  # elements remainder contains element numbers in interface
                raise ValueError(f'Number of lines ({elements0.shape[0]}) in {elem0_file} does not correspond with '
                                 f'the number of elements ({n_elem})')
            if self.timestep_start != 0:  # check if elements0 corresponds to timestep_start
                elem_file = join(self.dir_csm, f'CSM_Time{self.timestep_start}Surface{mp_id}Elements.dat')
                elements = np.loadtxt(elem_file)
                if int(elements[0, 0]) != n_elem or int(elements[0, 1]) != n_lp:
                    raise ValueError(f'Number of load points has changed for {mp_name}')

            # read in faces file for load points
            faces0_file = join(self.dir_csm, f'CSM_Time0Surface{mp_id}Cpu0Faces.dat')
            faces0 = np.loadtxt(faces0_file)
            if faces0.shape[1] != self.dimensions + 2:
                raise ValueError(f'Given dimension does not match coordinates')

            # get load point coordinates and ids of load points
            prev_elem = 0
            prev_lp = 0
            ids = np.arange(n_lp)
            coords_tmp = np.zeros((n_lp, 3))  # z-coordinate mandatory: 0.0 for 2D
            for i in range(0, n_lp):
                elem = int(faces0[i, 0])
                lp = int(faces0[i, 1])
                if elem < prev_elem:
                    raise ValueError(f'Element sequence is wrong ({elem}<{prev_elem})')
                elif elem == prev_elem and lp != prev_lp + 1:
                    raise ValueError(f'Next line for same element ({elem}) does not contain next load point')
                elif elem > prev_elem and lp != 1:
                    raise ValueError(f'First line for element ({elem}) does not contain its first load point')

                coords_tmp[i, :self.dimensions] = faces0[i, -self.dimensions:]  # extract last 'dimensions' columns

                prev_elem = elem
                prev_lp = lp

            x0 = coords_tmp[:, 0]
            y0 = coords_tmp[:, 1]
            z0 = coords_tmp[:, 2]

            # create ModelPart
            self.model.create_model_part(mp_name, x0, y0, z0, ids)

        # create output ModelParts (geometrical nodes)
        for item in (self.settings['interface_output']):
            mp_name = item['model_part']

            for i, surfaceID in enumerate(self.surfaceIDs):  # identify surfaceID corresponding to ModelPart
                if surfaceID in mp_name:
                    self.model_part_surface_ids[mp_name] = i
                    break
            if mp_name not in self.model_part_surface_ids:
                raise AttributeError(f'Could not identify surfaceID corresponding to ModelPart {mp_name}')
            mp_id = self.model_part_surface_ids[mp_name]

            # read in nodes file
            nodes0_file = join(self.dir_csm, f'CSM_Time0Surface{mp_id}Nodes.dat')
            nodes0 = np.loadtxt(nodes0_file, skiprows=1)   # first line is a header
            n_nodes0 = nodes0.shape[0]
            if nodes0.shape[1] != self.dimensions:
                raise ValueError(f'Given dimension does not match coordinates')
            if self.timestep_start != 0:  # check if nodes0 corresponds to timestep_start
                nodes_file = join(self.dir_csm, f'CSM_Time{self.timestep_start}Surface{mp_id}Nodes.dat')
                nodes = np.loadtxt(nodes_file, skiprows=1)  # first line is a header
                n_nodes = nodes.shape[0]
                if n_nodes != n_nodes0:
                    raise ValueError(f'Number of interface nodes has changed for {mp_name}')

            # get geometrical node coordinates
            ids = np.arange(n_nodes0)  # Abaqus does not use node ids but maintains the output order
            coords_tmp = np.zeros((n_nodes0, 3))  # z-coordinate mandatory: 0.0 for 2D
            coords_tmp[:, :self.dimensions] = nodes0

            x0 = coords_tmp[:, 0]
            y0 = coords_tmp[:, 1]
            z0 = coords_tmp[:, 2]

            # create ModelPart
            self.model.create_model_part(mp_name, x0, y0, z0, ids)

        # check whether the input ModelParts and output ModelParts have proper overlap
        for surfaceID in self.surfaceIDs:
            for item in self.settings['interface_input']:
                mp_name = item['model_part']
                if surfaceID in mp_name:
                    mp_in = self.model.get_model_part(mp_name)
                    break
            for item in self.settings['interface_output']:
                mp_name = item['model_part']
                if surfaceID in mp_name:
                    mp_out = self.model.get_model_part(mp_name)
                    break
            tools.check_bounding_box(mp_in, mp_out)

        # create Interfaces
        self.interface_input = data_structure.Interface(self.settings['interface_input'], self.model)
        self.interface_output = data_structure.Interface(self.settings['interface_output'], self.model)

        # run time
        self.run_time = 0.0

        # debug
        self.debug = False  # set on True to save copy of input and output files in every iteration

    def initialize(self):
        super().initialize()

    def initialize_solution_step(self):
        super().initialize_solution_step()

        self.iteration = 0
        self.timestep += 1

    @tools.time_solve_solution_step
    def solve_solution_step(self, interface_input):
        self.iteration += 1

        # store incoming loads
        self.interface_input.set_interface_data(interface_input.get_interface_data())

        # write loads (from interface data to a file that will be read by USR.f)
        self.write_loads()

        # copy input data for debugging
        if self.debug:
            for dct in self.interface_input.parameters:
                mp_name = dct['model_part']
                mp_id = self.model_part_surface_ids[mp_name]
                file_name1 = join(self.dir_csm, f'CSM_Time{self.timestep}Surface{mp_id}Cpu0Input.dat')
                file_name2 = join(self.dir_csm, f'CSM_Time{self.timestep}Surface{mp_id}Cpu0Input'
                                                f'_Iter{self.iteration}.dat')
                cmd = f'cp {file_name1} {file_name2}'
                os.system(cmd)

        # run Abaqus and check for (licensing) errors
        self.print_log(f'\n### Time step {self.timestep}, iteration {self.iteration} ###')
        bool_completed = False
        attempt = 0
        while not bool_completed and attempt < 10000:
            attempt += 1
            if attempt > 1:
                tools.print_info(f'Warning attempt {attempt - 1} to run Abaqus failed, new attempt in one minute',
                                 layout='warning')
                time.sleep(60)
                tools.print_info(f'Starting attempt {attempt}')
            if self.timestep == 1:
                cmd = f'abaqus job=CSM_Time{self.timestep} input=CSM_Time{self.timestep - 1} ' \
                    f'cpus={self.cores} output_precision=full interactive >> {self.logfile} 2>&1'
                subprocess.run(cmd, shell=True, cwd=self.dir_csm, executable='/bin/bash')
            else:
                cmd = f'abaqus job=CSM_Time{self.timestep} oldjob=CSM_Time{self.timestep - 1} input=CSM_Restart ' \
                    f'cpus={self.cores} output_precision=full interactive >> {self.logfile} 2>&1'
                subprocess.run(cmd, shell=True, cwd=self.dir_csm, executable='/bin/bash')

            # check log for completion and or errors
            subprocess.run(f'tail -n 10 {self.logfile} > Temp_log.coco', shell=True, cwd=self.dir_csm,
                           executable='/bin/bash')
            log_tmp = os.path.join(self.dir_csm, 'Temp_log.coco')
            bool_lic = True
            with open(log_tmp, 'r') as fp:
                for line in fp:
                    if any(x in line for x in ['Licensing error', 'license error',
                                               'Error checking out Abaqus license']):
                        bool_lic = False
            if not bool_lic:
                tools.print_info('Abaqus licensing error', layout='fail')
            elif 'COMPLETED' in line:  # check final line for completed
                bool_completed = True
            elif bool_lic:  # final line did not contain 'COMPLETED' but also no licensing error detected
                raise RuntimeError(f'Abaqus did not complete, unclassified error, see {self.logfile} for extra '
                                   f'information')

            # append additional information to log file
            cmd = f'tail -n 23 CSM_Time{self.timestep}.msg | head -n 15 | sed -e \'s/^[ \\t]*//\' ' \
                f'>> {self.logfile} 2>&1'
            subprocess.run(cmd, shell=True, cwd=self.dir_csm, executable='/bin/bash')

        # write Abaqus output
        cmd = f'abaqus ./GetOutput.exe CSM_Time{self.timestep} 1 >> {self.logfile} 2>&1'
        subprocess.run(cmd, shell=True, cwd=self.dir_csm, executable='/bin/bash')

        # read Abaqus output data
        for dct in self.interface_output.parameters:
            mp_name = dct['model_part']

            # read in output file for surface nodes
            mp_id = self.model_part_surface_ids[mp_name]
            file_name = join(self.dir_csm, f'CSM_Time{self.timestep}Surface{mp_id}Output.dat')
            data = np.loadtxt(file_name, skiprows=1)

            # copy output data for debugging
            if self.debug:
                file_name2 = join(self.dir_csm, f'CSM_Time{self.timestep}Surface{mp_id}Output_Iter{self.iteration}.dat')
                shutil.copy(file_name, file_name2)

            if data.shape[1] != self.dimensions:
                raise ValueError(f'given dimension does not match coordinates')

            # get surface node displacements
            n_nodes = data.shape[0]
            model_part = self.model.get_model_part(mp_name)
            if n_nodes != model_part.size:
                raise ValueError('size of data does not match size of ModelPart')

            displacement = np.zeros((n_nodes, 3))  # also require z-input for 2D cases
            displacement[:, :self.dimensions] = data

            self.interface_output.set_variable_data(mp_name, 'displacement', displacement)

        return self.interface_output

    def finalize_solution_step(self):
        super().finalize_solution_step()
        if self.timestep and (self.timestep - 1) % self.settings['save_iterations']:
            to_be_removed_suffix = ['.com', '.dat', '.mdl', '.msg', '.odb', '.prt', '.res', '.sim', '.sta', '.stt',
                                    'Surface*Cpu0Input.dat', 'Surface*Output.dat']
            cmd = ''
            for suffix in to_be_removed_suffix:
                cmd += f'rm CSM_Time{self.timestep - 1}{suffix}; '
            subprocess.run(cmd, shell=True, cwd=self.dir_csm, executable='/bin/bash')

    def finalize(self):
        super().finalize()

    def get_interface_input(self):
        return self.interface_input

    def get_interface_output(self):
        return self.interface_output

    def check_software(self):
        """Check whether the software requirements for this wrapper are fulfilled."""
        # Python version: 3.6 or higher
        if sys.version_info < (3, 6):
            raise RuntimeError('Python version 3.6 or higher required')

        # Abaqus availability
        if shutil.which('abaqus') is None:
            raise RuntimeError('Abaqus must be available')

        # Abaqus version
        result = subprocess.run(['abaqus', 'information=release'], stdout=subprocess.PIPE)
        if self.version not in str(result.stdout):
            raise RuntimeError(f'Abaqus version {self.version} is required')

        # compilers
        if shutil.which('ifort') is None:
            raise RuntimeError('Intel compiler ifort must be available')

    def print_log(self, msg):
        with open(os.path.join(self.dir_csm, self.logfile), 'a') as f:
            print(msg, file=f)

    # noinspection PyMethodMayBeStatic
    def make_elements(self, face_file, output_file):
        """
        Reads the CSM_Time0SurfaceXFaces.dat file created by USRInit.f and converts it to
        CSM_Time0SurfaceXElements.dat, formatted such that USR.f can read it. The first line mentions the number of
        elements and load integration points. All lines below contain the element number and cumulative number of
        load integration points. The USR.f uses this to quickly find the correct line to read from the load input files.
        """
        face_file = Path(face_file)
        output_file = Path(output_file)
        first_loop = True
        element_prev = -1
        point_prev = -1
        element_0 = -1
        point_0 = -1
        n_elem = 0  # total number of elements
        n_lp = 0    # total number of load points
        element_list = []

        with open(face_file, 'r') as file:
            for num, line in enumerate(file, start=1):
                values = line.strip().split()
                element = int(values[0])
                point = int(values[1])
                if element == element_0 and point == point_0:  # all values multiple times in file, only once needed
                    break
                if element == element_prev:  # load points of same element
                    if point == point_prev + 1:  # deviation from ascending order indicates mistake
                        point_prev = point
                        n_lp += 1
                    else:
                        if point > point_prev:
                            msg = textwrap.fill(f'Error while processing {face_file.name} line {num}. Load point number'
                                                f' increases by more 1 one per line for element {element}, found '
                                                f'{point} but {point_prev + 1} was expected', width=80)
                            raise ValueError(msg)
                        elif point == 1:
                            msg = textwrap.fill(f'Error while processing {face_file.name} line {num}. Load point number'
                                                f' {point} is lower than previous ({point_prev}). Check if surface '
                                                f'contains element with multiple faces in surface, for example at '
                                                f'corners. This is not allowed.', width=80)
                            raise ValueError(msg)
                        else:
                            msg = textwrap.fill(f'Error while processing {face_file.name} line {num}. Load point number'
                                                f' {point} is lower than previous ({point_prev}).', width=80)
                            raise ValueError(msg)
                else:  # new element started
                    if point == 1:
                        # add to elements file the element number
                        # and cumulative number of load points encountered BEFORE reaching this element
                        element_list.append([element, n_lp])
                        point_prev = point
                        element_prev = element
                        n_elem += 1
                        n_lp += 1
                        if first_loop:
                            element_0 = element
                            point_0 = point
                            first_loop = False
                    else:
                        msg = textwrap.fill(
                            f'Error while processing {face_file.name} line {num}: Load point number does not start at 1'
                            f' for element {element}, {point} was found instead.', width=80)
                        raise ValueError(msg)

        element_list.insert(0, [n_elem, n_lp])
        element_list = np.array(element_list)
        np.savetxt(output_file, element_list, fmt='%10d')

    # noinspection PyMethodMayBeStatic
    def replace_fortran(self, line, orig, new):
        """
        The length of a line in FORTRAN 77 (i.e. fixed-form) is limited, replacing working directories can exceed this
        limit. This functions splits these strings over multiple lines.
        """
        if '|' in line:
            char_limit = 72
            wrapper = textwrap.TextWrapper(width=char_limit, subsequent_indent='     &', replace_whitespace=False,
                                           drop_whitespace=False)
            line = wrapper.fill(line.replace(orig, new))
        return line

    def write_start_and_restart_inp(self, input_file, output_file, restart_file):
        """
        Read the case-file (.inp) and process it in an input-file for the first time step and a restart file used for
        all subsequent time steps.
        """
        bool_restart = False

        rf = open(restart_file, 'w')
        of = open(output_file, 'w')

        rf.write('*HEADING \n')
        rf.write('*RESTART, READ \n')

        with open(input_file) as f:
            line = f.readline()
            while line:
                if '*step' in line.lower():
                    contents = line.split(',')  # split string on commas
                    line_new = ''
                    for s in contents:
                        if s.strip().startswith('inc='):
                            numbers = re.findall(r'\d+', s)
                            if (not self.subcycling) and int(numbers[0]) != 1:
                                raise NotImplementedError(f'inc={numbers[0]}: subcycling was not requested but '
                                                          f'max_num_inc > 1')
                            else:
                                line_new += f' inc={self.max_num_inc},'
                        else:
                            line_new += s+','
                    line_new = line_new[:-1]+'\n'  # remove the last added comma and add a newline

                    of.write(line_new)
                    if bool_restart:
                        rf.write(line_new)
                    line = f.readline()
                elif '*dynamic' in line.lower():
                    contents = line.split(',')
                    for s in contents:
                        if 'application' in s.lower():
                            contents_B = s.split('=')
                            app = contents_B[1].lower().strip()
                            if app == 'quasi-static' or app == 'moderate dissipation' or 'transient fidelity':
                                if not self.subcycling:
                                    line_2 = f'{self.delta_t}, {self.delta_t},\n'
                                else:
                                    line_2 = f'{self.initial_inc}, {self.delta_t}, {self.min_inc}, {self.max_inc}\n'
                            else:
                                raise NotImplementedError(
                                    f'{contents_B[1]} not available: currently only quasi-static, moderate dissipation '
                                    f'and transient fidelity are implemented for the Abaqus wrapper')
                    of.write(line)
                    if bool_restart:
                        rf.write(line)
                        f.readline()  # need to skip the next line
                        of.write(line_2)  # change the time step in the Abaqus step
                        if bool_restart:
                            rf.write(line_2)  # change the time step in the Abaqus step (restart-file)
                        line = f.readline()
                    elif '*static' in line.lower():
                        of.write(line)
                        if bool_restart:
                            rf.write(line)
                        f.readline()  # need to skip the next line
                        if not self.subcycling:
                            raise NotImplementedError(
                                'Only Static with subcycling is implemented for the Abaqus wrapper')
                        line_2 = f'{self.initial_inc}, {self.delta_t}, {self.min_inc}, {self.max_inc}\n'
                        f.readline()
                    of.write(line_2)  # change the time step in the Abaqus step
                    if bool_restart:
                        rf.write(line_2)  # change the time step in the Abaqus step (restart-file)
                    line = f.readline()
                else:
                    of.write(line)
                    if bool_restart:
                        rf.write(line)
                    line = f.readline()
                if '** --'in line:
                    bool_restart = True
        rf.close()
        of.close()

    def write_loads(self):
        """Write the incoming loads to a file that is read by the READDATA subroutine in USR.f."""
        for dct in self.interface_input.parameters:
            mp_name = dct['model_part']
            mp_id = self.model_part_surface_ids[mp_name]
            model_part = self.model.get_model_part(mp_name)

            pressure = self.interface_input.get_variable_data(mp_name, 'pressure')
            traction = self.interface_input.get_variable_data(mp_name, 'traction')
            data = np.hstack((pressure, traction[:, :self.dimensions]))
            fmt = (self.dimensions + 1) * '%27.17e'  # format of load input file should be exactly this for USR.f
            file_name = join(self.dir_csm, f'CSM_Time{self.timestep}Surface{mp_id}Cpu0Input.dat')
            np.savetxt(file_name, data, fmt=fmt, header=f'{model_part.size}', comments='')

            # start of a simulation with ramp, needs an initial load at time 0: set at zero load
            if self.iteration == 1 and self.timestep == 1 and self.ramp:
                file_name = join(self.dir_csm, f'CSM_Time0Surface{mp_id}Cpu0Input.dat')
                np.savetxt(file_name, data * 0, fmt=fmt, header=f'{model_part.size}', comments='')
