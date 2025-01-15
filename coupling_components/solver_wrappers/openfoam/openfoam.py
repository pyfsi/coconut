from coconut import data_structure
from coconut.coupling_components.solver_wrappers.solver_wrapper import SolverWrapper
from coconut.data_structure.interface import Interface
from coconut import tools
from coconut.coupling_components.solver_wrappers.openfoam import openfoam_io as of_io

import numpy as np
import os
import shutil
import time
import subprocess
import re
from glob import glob
from subprocess import check_call


def create(parameters):
    return SolverWrapperOpenFOAM(parameters)


class SolverWrapperOpenFOAM(SolverWrapper):
    version = None  # OpenFOAM version with dot, e.g. 8 , set in subclass
    check_coupling_convergence_possible = True  # can solver check convergence after 1 iteration?

    # define input and output variables
    accepted_in_var = ['displacement']
    accepted_out_var = ['pressure', 'traction']

    @tools.time_initialize
    def __init__(self, parameters):
        super().__init__(parameters)

        if self.version is None:
            raise NotImplementedError(
                'Base class method called, class variable version needs to be set in the derived class')

        # settings
        self.settings = parameters['settings']
        self.working_directory = self.settings['working_directory']
        self.env = None  # environment in which correct version of OpenFOAM software is available, set in subclass

        # adapted application from openfoam ('coconut_<application name>')
        self.application = self.settings['application']
        self.delta_t = self.settings['delta_t']
        self.time_precision = self.settings['time_precision']
        self.timestep_start = self.settings['timestep_start']
        self.start_time = self.timestep_start * self.delta_t
        self.timestep = self.physical_time = self.iteration = self.prev_timestamp = self.cur_timestamp = None
        self.save_restart = self.settings['save_restart']
        self.openfoam_process = None
        self.write_interval = self.write_precision = None
        self.fext = None
        self.nheaderfooter = None

        # boundary_names is the set of boundaries in OpenFoam used for coupling
        self.boundary_names = self.settings['boundary_names']
        self.cores = None
        self.model = None

        # set on True if you want to clean the adapted application and compile.
        self.compile_clean = self.settings.get('compile_clean', False)

        # remove possible CoCoNuT-message from previous interrupt
        self.coco_messages = tools.CocoMessages(self.working_directory, max_wait_time=600,
                                                timed_out_action=self.timed_out)
        self.coco_messages.remove_all_messages()

        # variables for kinematic pressure and traction conversion
        self.density_for_pressure = None
        self.density_for_traction = None
        self.wall_shear_stress_variable = None
        self.wall_shear_stress_function_object_library = None

        # residual variables
        self.residual_variables = self.settings.get('residual_variables', None)
        self.res_filepath = os.path.join(self.working_directory, 'residuals.csv')
        self.mp_in_decompose_seq_dict = {}
        self.mp_out_reconstruct_seq_dict = {}

        if self.residual_variables is not None:
            self.write_residuals_fileheader()

    @tools.time_initialize
    def initialize(self):
        super().initialize()

        # check interface names in 'interface_input' and 'interface_output' with boundary names provided in
        # boundary_names
        self.check_interfaces()

        # obtain number of cores from self.working_directory/system/decomposeParDict
        self.cores = 1
        if self.settings['parallel']:
            file_name = os.path.join(self.working_directory, 'system/decomposeParDict')
            if not os.path.isfile(file_name):
                raise RuntimeError(
                    f'In the JSON parameters "parallel" is set to {True} but {file_name} does not exist')
            else:
                with open(file_name, 'r') as file:
                    decomposedict_string = file.read()
                self.cores = of_io.get_int(input_string=decomposedict_string, keyword='numberOfSubdomains')

        # take density into account for solvers working with kinematic pressure or traction
        self.kinematic_conversion()

        # modify controlDict file to add pressure and wall shear stress functionObjects for all the boundaries in
        # self.settings["boundary_names"]
        self.read_modify_controldict()

        if self.write_interval is not None and self.save_restart % self.write_interval:
            raise RuntimeError(
                f'self.save_restart (= {self.save_restart}) should be an integer multiple of writeInterval '
                f'(= {self.write_interval}). Modify the controlDict accordingly.')

        # creating Model
        self.model = data_structure.Model()

        # writeCellcentres writes cellcentres in internal field and face centres in boundaryField
        self.write_cell_centres()

        boundary_filename = os.path.join(self.working_directory, 'constant/polyMesh/boundary')
        for boundary in self.boundary_names:
            with open(boundary_filename, 'r') as boundary_file:
                boundary_file_string = boundary_file.read()
            boundary_dict = of_io.get_dict(input_string=boundary_file_string, keyword=boundary)
            # get point ids and coordinates for all the faces in the boundary
            node_ids, node_coords = of_io.get_boundary_points(case_directory=self.working_directory, time_folder='0',
                                                              boundary_name=boundary)
            nfaces = of_io.get_int(input_string=boundary_dict, keyword='nFaces')
            start_face = of_io.get_int(input_string=boundary_dict, keyword='startFace')

            # create input model part
            self.model.create_model_part(f'{boundary}_input', node_coords[:, 0], node_coords[:, 1], node_coords[:, 2],
                                         node_ids)

            x0, y0, z0 = self.read_face_centres(boundary, nfaces)
            ids = np.arange(start_face, start_face + nfaces)

            # create output model part
            self.model.create_model_part(f'{boundary}_output', x0, y0, z0, ids)

        # create interfaces
        self.interface_input = Interface(self.settings['interface_input'], self.model)
        self.interface_output = Interface(self.settings['interface_output'], self.model)

        # define timestep and physical time
        self.timestep = self.timestep_start
        self.physical_time = self.start_time

        # if parallel do a decomposition and establish a remapping for the output based on the faceProcAddressing
        # Note concerning the sequence: The file ./processorX/constant/polyMesh/faceprocAddressing contains a list of
        # indices referring to the original index in the ./constant/polyMesh/faces file, these indices go from 0 to
        # nfaces -1
        # However, mesh faces can be shared between processors and it has to be tracked whether these are inverted
        # or not
        # This inversion is indicated by negative indices
        # However, as minus 0 is not a thing, the indices are first incremented by 1 before inversion
        # Therefore to get the correct index one should use |index|-1!!

        if self.settings['parallel']:
            if self.start_time == 0:
                subprocess.check_call(f'decomposePar -force -time {self.start_time} &> log.decomposePar',
                                      cwd=self.working_directory, shell=True, env=self.env)

            for boundary in self.boundary_names:
                mp_out_name = f'{boundary}_output'
                mp_output = self.model.get_model_part(f'{boundary}_output')
                nfaces = mp_output.size
                start_face = mp_output.id[0]
                self.mp_out_reconstruct_seq_dict[mp_out_name] = []
                for p in range(self.cores):
                    path = os.path.join(self.working_directory, f'processor{p}/constant/polyMesh/faceProcAddressing')
                    with open(path, 'r') as f:
                        face_proc_add_string = f.read()
                    face_proc_add = np.abs(of_io.get_scalar_array(input_string=face_proc_add_string, is_int=True))
                    face_proc_add -= 1  # in openfoam face ids are incremented by 1
                    self.mp_out_reconstruct_seq_dict[mp_out_name] += (face_proc_add[(face_proc_add >= start_face) & (
                            face_proc_add < start_face + nfaces)] - start_face).tolist()

                if len(self.mp_out_reconstruct_seq_dict[mp_out_name]) != nfaces:
                    print(f'sequence: {len(mp_output.sequence)}')
                    print(f'nNodes: {mp_output.size}')
                    raise ValueError('Number of face indices in sequence does not correspond to number of faces')

                mp_in_name = f'{boundary}_input'
                mp_input = self.model.get_model_part(mp_in_name)
                self.mp_in_decompose_seq_dict[mp_in_name] = {}
                # get the point sequence in the boundary for points in different processors
                for p in range(self.cores):
                    proc_dir = os.path.join(self.working_directory, f'processor{p}')
                    point_ids, points = of_io.get_boundary_points(proc_dir, '0', boundary)

                    if point_ids.size:
                        with open(os.path.join(proc_dir, 'constant/polyMesh/pointProcAddressing'), 'r') as f:
                            point_proc_add = np.abs(of_io.get_scalar_array(input_string=f.read(), is_int=True))
                        sorter = np.argsort(mp_input.id)
                        self.mp_in_decompose_seq_dict[mp_in_name][p] = sorter[
                            np.searchsorted(mp_input.id, point_proc_add[point_ids], sorter=sorter)]
                    else:
                        self.mp_in_decompose_seq_dict[mp_in_name][p] = None

        # starting the OpenFOAM infinite loop for coupling!
        if not self.settings['parallel']:
            cmd = self.application + '&> log.' + self.application
        else:
            cmd = 'mpirun -np ' + str(self.cores) + ' ' + self.application + ' -parallel &> log.' + self.application

        self.openfoam_process = subprocess.Popen(cmd, cwd=self.working_directory, shell=True, env=self.env)

        # pass on process to coco_messages for polling
        self.coco_messages.set_process(self.openfoam_process)

    def initialize_solution_step(self):
        super().initialize_solution_step()

        timestamp = '{:.{}f}'.format(self.physical_time, self.time_precision)

        # prepare new time step folder and reset the number of iterations
        self.timestep += 1
        self.iteration = 0
        self.physical_time += self.delta_t

        self.prev_timestamp = timestamp
        self.cur_timestamp = f'{self.physical_time:.{self.time_precision}f}'

        if not self.settings['parallel']:  # if serial
            new_path = os.path.join(self.working_directory, self.cur_timestamp)
            if os.path.isdir(new_path):
                tools.print_info(f'Overwrite existing time step folder: {new_path}', layout='warning')
                subprocess.check_call(f'rm -rf {new_path}', shell=True)
        else:
            for i in np.arange(self.cores):
                new_path = os.path.join(self.working_directory, 'processor' + str(i), self.cur_timestamp)
                if os.path.isdir(new_path):
                    if i == 0:
                        tools.print_info(f'Overwrite existing time step folder: {new_path}', layout='warning')
                    subprocess.check_call(f'rm -rf {new_path}', shell=True)

        self.coco_messages.send_message('next')
        self.coco_messages.wait_message('next_ready')

    @tools.time_solve_solution_step
    def solve_solution_step(self, interface_input):
        self.iteration += 1

        # store incoming displacements
        self.interface_input.set_interface_data(interface_input.get_interface_data())

        # write interface data to OpenFOAM-file
        self.write_node_input()

        # copy output data for debugging
        if self.debug:
            if self.cores > 1:
                for i in range(0, self.cores):
                    path_from = os.path.join(self.working_directory, f'processor{i}/constant/pointDisplacementTmp')
                    path_to = os.path.join(self.working_directory,
                                           f'processor{i}/constant/'
                                           f'pointDisplacementTmp_{self.timestep}_{self.iteration}')
                    shutil.copy(path_from, path_to)
            else:
                path_from = os.path.join(self.working_directory, 'constant/pointDisplacementTmp')
                path_to = os.path.join(self.working_directory,
                                       f'constant/pointDisplacementTmp_{self.timestep}_{self.iteration}')
                shutil.copy(path_from, path_to)

        self.delete_prev_iter_output()

        self.coco_messages.send_message('continue')
        self.coco_messages.wait_message('continue_ready')

        if self.check_coupling_convergence:
            # check if OpenFOAM converged after 1 iteration
            self.coco_messages.wait_message('check_ready')
            self.coupling_convergence = self.coco_messages.check_message('solver_converged')
            if self.print_coupling_convergence and self.coupling_convergence:
                tools.print_info(f'{self.__class__.__name__} converged')

        # read data from OpenFOAM
        self.read_node_output()

        # copy output data for debugging
        if self.debug:
            for boundary in self.boundary_names:
                post_process_time_folder = os.path.join(self.working_directory,
                                                        f'postProcessing/coconut_{boundary}/surface',
                                                        self.cur_timestamp)
                for file in (f'{self.wall_shear_stress_variable}_patch_{boundary}{self.fext}',
                             f'p_patch_{boundary}{self.fext}'):
                    path = os.path.join(post_process_time_folder, file)
                    iter_path = os.path.join(post_process_time_folder,
                                             f'{file[:-len(self.fext)]}_{self.iteration}{self.fext}')
                    shutil.copy(path, iter_path)

        # return interface_output object
        return self.interface_output

    def finalize_solution_step(self):
        super().finalize_solution_step()

    @tools.time_save
    def output_solution_step(self):
        super().output_solution_step()

        if not self.debug:
            for boundary in self.boundary_names:
                post_process_time_folder = os.path.join(self.working_directory,
                                                        f'postProcessing/coconut_{boundary}/surface',
                                                        self.cur_timestamp)
                shutil.rmtree(post_process_time_folder)

        self.coco_messages.send_message('save')
        self.coco_messages.wait_message('save_ready')

        if self.residual_variables is not None:
            self.write_of_residuals()

    def finalize(self):
        super().finalize()

        self.coco_messages.send_message('stop')
        self.coco_messages.wait_message('stop_ready')
        self.openfoam_process.wait()
        if not self.debug:
            files_to_delete = glob(os.path.join(self.working_directory, 'constant/pointDisplacementTmp*')) + glob(
                os.path.join(self.working_directory, 'processor*/constant/pointDisplacementTmp*'))
            for filepath in files_to_delete:
                os.remove(filepath)

            for boundary in self.boundary_names:
                post_process_folder = os.path.join(self.working_directory, f'postProcessing/coconut_{boundary}')
                shutil.rmtree(post_process_folder, ignore_errors=True)

    def compile_adapted_openfoam_solver(self):
        # compile openfoam adapted solver
        solver_dir = os.path.join(os.path.dirname(__file__), f'v{self.version.replace(".", "")}', self.application)
        try:
            if self.compile_clean:
                subprocess.check_call(f'wclean {solver_dir} && wmake {solver_dir} &> log.wmake',
                                      cwd=self.working_directory, shell=True,
                                      env=self.env)
            else:
                subprocess.check_call(f'wmake {solver_dir} &> log.wmake', cwd=self.working_directory, shell=True,
                                      env=self.env)
        except subprocess.CalledProcessError:
            raise RuntimeError(
                f'Compilation of {self.application} failed. Check {os.path.join(self.working_directory, "log.wmake")}')

    def write_cell_centres(self):
        check_call('postProcess -func writeCellCentres -time 0 &> log.writeCellCentres;', cwd=self.working_directory,
                   shell=True, env=self.env)

    def read_face_centres(self, boundary_name, nfaces):
        raise NotImplementedError('Base class method is called, should be implemented in derived class')

    def delete_prev_iter_output(self):
        # pressure and wall shear stress files are removed to avoid openfoam to append data in the new iteration
        for boundary in self.boundary_names:
            # specify location of pressure and traction
            post_process_time_folder = os.path.join(self.working_directory,
                                                    f'postProcessing/coconut_{boundary}/surface', self.cur_timestamp)
            for file in (f'{self.wall_shear_stress_variable}_patch_{boundary}{self.fext}',
                         f'p_patch_{boundary}{self.fext}'):
                path = os.path.join(post_process_time_folder, file)
                if os.path.isfile(path):
                    os.remove(path)

    def read_node_output(self):
        """
        reads the pressure and wall shear stress from the <case directory>/postProcessing for serial and parallel. In
        case of parallel, it uses mp.sequence (using faceProcAddressing) to map the values to the face centres.

        :return:
        """
        for boundary in self.boundary_names:
            # specify location of pressure and traction
            mp_name = f'{boundary}_output'
            mp = self.model.get_model_part(mp_name)
            nfaces = mp.size
            post_process_time_folder = os.path.join(self.working_directory,
                                                    f'postProcessing/coconut_{boundary}/surface', self.cur_timestamp)
            wss_filename = os.path.join(post_process_time_folder,
                                        f'{self.wall_shear_stress_variable}_patch_{boundary}{self.fext}')
            pres_filepath = os.path.join(post_process_time_folder, f'p_patch_{boundary}{self.fext}')

            # check if the pressure and wall shear stress files completed by openfoam and read data
            self.check_output_file(wss_filename, nfaces)
            wss_tmp = np.loadtxt(wss_filename, comments='#')[:, 3:]
            self.check_output_file(pres_filepath, nfaces)
            pres_tmp = np.loadtxt(pres_filepath, comments='#')[:, 3]

            if self.settings['parallel']:
                pos_list = self.mp_out_reconstruct_seq_dict[mp_name]
            else:
                pos_list = [pos for pos in range(0, nfaces)]

            wall_shear_stress = np.empty_like(wss_tmp)
            pressure = np.empty((pres_tmp.size, 1))

            wall_shear_stress[pos_list] = wss_tmp[:, ]
            pressure[pos_list, 0] = pres_tmp

            self.interface_output.set_variable_data(mp_name, 'traction', -wall_shear_stress * self.density_for_traction)
            self.interface_output.set_variable_data(mp_name, 'pressure', pressure * self.density_for_pressure)

    # noinspection PyMethodMayBeStatic
    def write_footer(self, file_name):
        # write OpenFOAM-footer at the end of file
        with open(file_name, 'a') as f:
            f.write('\n// ************************************************************************* //\n')

    def write_node_input(self):
        """
        creates pointDisplacementTmp for supplying the displacement field in the FSI coupling. This file is created
        from 0/pointDisplacement. The boundary field for boundaries participating in the FSI coupling is modified to
        supply the boundary displacement field from structural solver. If the OpenFOAM solver is run in parallel,
        the field is subsequently decomposed using the command: decomposePar.
       :return:
       """
        if not self.settings['parallel']:
            pointdisp_filename_ref = os.path.join(self.working_directory, '0/pointDisplacement')
            pointdisp_filename = os.path.join(self.working_directory, 'constant/pointDisplacementTmp')

            with open(pointdisp_filename_ref, 'r') as ref_file:
                pointdisp_string = ref_file.read()

            for boundary in self.boundary_names:
                mp_name = f'{boundary}_input'
                displacement = self.interface_input.get_variable_data(mp_name, 'displacement')
                boundary_dict = of_io.get_dict(input_string=pointdisp_string, keyword=boundary)
                boundary_dict_new = of_io.update_vector_array_dict(dict_string=boundary_dict, vector_array=displacement)
                pointdisp_string = pointdisp_string.replace(boundary_dict, boundary_dict_new)

            with open(pointdisp_filename, 'w') as f:
                f.write(pointdisp_string)
        else:
            for proc in range(self.cores):
                pointdisp_filename_ref = os.path.join(self.working_directory, f'processor{proc}', '0/pointDisplacement')
                pointdisp_filename = os.path.join(self.working_directory, f'processor{proc}',
                                                  'constant/pointDisplacementTmp')

                with open(pointdisp_filename_ref, 'r') as ref_file:
                    pointdisp_string = ref_file.read()

                for boundary in self.boundary_names:
                    mp_name = f'{boundary}_input'
                    displacement = self.interface_input.get_variable_data(mp_name, 'displacement')
                    seq = self.mp_in_decompose_seq_dict[mp_name][proc]
                    if seq is not None:
                        boundary_dict = of_io.get_dict(input_string=pointdisp_string, keyword=boundary)
                        boundary_dict_new = of_io.update_vector_array_dict(dict_string=boundary_dict,
                                                                           vector_array=displacement[seq])
                        pointdisp_string = pointdisp_string.replace(boundary_dict, boundary_dict_new)

                with open(pointdisp_filename, 'w') as f:
                    f.write(pointdisp_string)

    # old way of implementation (with decomposePar)
    # if self.settings['parallel']:
    #     subprocess.check_call(f'decomposePar -fields -time '' -constant &> log.decomposePar;',
    #                cwd=self.working_directory, shell=True, env=self.env)

    # noinspection PyMethodMayBeStatic
    def check_output_file(self, filename, nfaces):
        counter = 0
        nlines = 0
        lim = 10000
        sleep_time = 0.01
        while (nlines < nfaces + self.nheaderfooter) and counter < lim:
            if os.path.isfile(filename):
                with open(filename, 'r') as f:
                    nlines = sum(1 for _ in f)
            time.sleep(sleep_time)
            counter += 1
            if not counter % 1000:
                tools.print_info(f'Waiting {counter * sleep_time} s for {filename}')
        if counter == lim:
            raise RuntimeError(f'Timed out waiting for file: {filename}')
        else:
            return True

    def timed_out(self, message):
        self.openfoam_process.kill()
        self.openfoam_process.wait()
        raise RuntimeError(f'CoCoNuT timed out in the OpenFOAM solver_wrapper, waiting for message: {message}.coco')

    def check_software(self):
        try:
            version_nr = subprocess.check_output('echo $WM_PROJECT_VERSION',
                                                 shell=True, env=self.env).decode('utf-8').strip()  # output b'XX\n'
        except subprocess.CalledProcessError:
            raise RuntimeError(
                f'OpenFOAM not loaded properly. Check if the solver load commands for the "machine_name"'
                f' are correct in solver_modules.py.')

        # check version
        if version_nr != self.version:
            raise RuntimeError(
                f'OpenFOAM-{self.version} should be loaded! Currently, OpenFOAM-{version_nr} is loaded.'
                f' Check if the solver load commands for the "machine_name" are correct in solver_modules.py.')

    def check_interfaces(self):
        """
        checks the dictionaries from 'interface_input' and 'interface_output' in parameters.json file. The model part
        name must be the concatenation of an entry from `boundary_names` and the string `_input`, for 'interface_input'
        and for 'interface_output' it must be the concatenation of an entry from `boundary_names` and the string
        `_output`.
        :return:
        """
        input_interface_model_parts = [param['model_part'] for param in self.settings['interface_input']]
        output_interface_model_parts = [param['model_part'] for param in self.settings['interface_output']]
        boundary_names = self.settings['boundary_names']

        for boundary_name in boundary_names:
            if f'{boundary_name}_input' not in input_interface_model_parts:
                raise RuntimeError(
                    f'Error in json file: {boundary_name}_input not listed in "interface_input": '
                    f'{self.settings["interface_input"]}.\n. <boundary> in the "boundary_names" in json file should '
                    f'have corresponding <boundary>_input in "interface_input" list')

            if f'{boundary_name}_output' not in output_interface_model_parts:
                raise RuntimeError(
                    f'Error in json file: {boundary_name}_output not listed in "interface_output": '
                    f'{self.settings["interface_output"]}.\n. <boundary> in the "boundary_names" in json file should '
                    f'have corresponding <boundary>_output in "interface_output" list')

    def read_modify_controldict(self):
        """
        reads the controlDict file in the case-directory and modifies some entries required by the coconut_pimpleFoam.
        The values of these entries are taken from paramters.json file.
        :return:
        """

        file_name = os.path.join(self.working_directory, 'system/controlDict')
        with open(file_name, 'r') as control_dict_file:
            control_dict = control_dict_file.read()
        write_control = of_io.get_string(input_string=control_dict, keyword='writeControl')
        if write_control == 'timeStep':
            self.write_interval = of_io.get_int(input_string=control_dict, keyword='writeInterval')
        time_format = of_io.get_string(input_string=control_dict, keyword='timeFormat')
        self.write_precision = of_io.get_int(input_string=control_dict, keyword='writePrecision')

        if not time_format == 'fixed':
            msg = f'timeFormat:{time_format} in controlDict not implemented. Changed to "fixed"'
            tools.print_info(msg, layout='warning')
            control_dict = re.sub(r'timeFormat' + of_io.delimiter + r'\w+', f'timeFormat    fixed',
                                  control_dict)
        control_dict = re.sub(r'application' + of_io.delimiter + r'\w+', f'{"application":<16}{self.application}',
                              control_dict)
        control_dict = re.sub(r'startTime' + of_io.delimiter + of_io.float_pattern,
                              f'{"startTime":<16}{self.start_time}', control_dict)
        control_dict = re.sub(r'deltaT' + of_io.delimiter + of_io.float_pattern, f'{"deltaT":<16}{self.delta_t}',
                              control_dict)
        control_dict = re.sub(r'timePrecision' + of_io.delimiter + of_io.int_pattern,
                              f'{"timePrecision":<16}{self.time_precision}',
                              control_dict)
        control_dict = re.sub(r'endTime' + of_io.delimiter + of_io.float_pattern, f'{"endTime":<16}1e15', control_dict)

        # delete previously defined coconut functions
        coconut_start_string = '\n\n// CoCoNuT function objects'
        control_dict = re.sub(coconut_start_string + r'.*', '', control_dict, flags=re.S)

        with open(file_name, 'w') as control_dict_file:
            control_dict_file.write(control_dict)
            control_dict_file.write(coconut_start_string + '\n')

            control_dict_file.write('boundaryNames   (' + ' '.join(self.boundary_names) + ');\n\n')

            if self.check_coupling_convergence:
                control_dict_file.write(f'checkCouplingConvergence    true;\n\n')

            control_dict_file.write('functions\n{\n')
            dct, name = self.wall_shear_stress_dict()
            control_dict_file.write(dct)
            function_object_names = [name]
            for boundary_name in self.boundary_names:
                dct, name = self.pressure_and_traction_dict(boundary_name)
                control_dict_file.write(dct)
                function_object_names.append(name)
            control_dict_file.write('}\n\n')

            control_dict_file.write('coconutFunctionObjects \n(\n\t' + '\n\t'.join(function_object_names)+'\n);\n')

            self.write_footer(file_name)

    def wall_shear_stress_dict(self):
        raise NotImplementedError('Base class method is called, should be implemented in derived class')

    def pressure_and_traction_dict(self, boundary_name):
        raise NotImplementedError('Base class method is called, should be implemented in derived class')

    def kinematic_conversion(self):
        raise NotImplementedError('Base class method is called, should be implemented in derived class')

    def write_residuals_fileheader(self):
        header = ''
        sep = ', '
        with open(self.res_filepath, 'w') as f:
            f.write('# Residuals\n')
            for variable in self.residual_variables:
                header += variable + sep
            f.write(header.strip(sep) + '\n')

    def write_of_residuals(self):
        """
        it reads the log file generated by coconut_pimpleFoam solver and writes the last initial residual of the fields
        in the pimple iterations, for every coupling iteration. The fields should be given in the parameters.json file
        with the key-'residual_variables' and values-list(OpenFOAM variables), e.g. 'residual_variables': ['U', 'p']
        :return:
        """
        log_filepath = os.path.join(self.working_directory, f'log.{self.application}')
        if os.path.isfile(log_filepath):
            with open(log_filepath, 'r') as f:
                log_string = f.read()
            time_start_string = f'Time = {self.prev_timestamp}'
            time_end_string = f'Time = {self.cur_timestamp}'
            match = re.search(time_start_string + r'(.*)' + time_end_string, log_string, flags=re.S)
            if match is not None:
                time_block = match.group(1)
                iteration_block_list = re.findall(
                    r'Coupling iteration = \d+(.*?)Coupling iteration \d+ end', time_block, flags=re.S)
                for iteration_block in iteration_block_list:
                    residual_array = np.empty(len(self.residual_variables))
                    for i, variable in enumerate(self.residual_variables):
                        search_string = f'Solving for {variable}, Initial residual = ({of_io.float_pattern})'
                        var_residual_list = re.findall(search_string, iteration_block)
                        if var_residual_list:
                            # last initial residual of pimple loop
                            var_residual = float(var_residual_list[-1])
                            residual_array[i] = var_residual
                        else:
                            raise RuntimeError(f'Variable: {variable} equation is not solved')

                    with open(self.res_filepath, 'a') as f:
                        np.savetxt(f, [residual_array], delimiter=', ')
