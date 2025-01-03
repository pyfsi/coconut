import glob
import os
import re
import shutil
import subprocess
import time
import warnings
import xml.etree.ElementTree as ElTr
from getpass import getuser
from os.path import join
from xml.dom import minidom

import numpy as np
from coconut import data_structure
from coconut import tools
from coconut.coupling_components.solver_wrappers.solver_wrapper import SolverWrapper


def create(parameters):
    return SolverWrapperAbaqusCSE(parameters)


class SolverWrapperAbaqusCSE(SolverWrapper):
    version = None  # Abaqus version, e.g. 2024, set in subclass
    check_coupling_convergence_possible = False  # can solver check convergence after 1 iteration?

    # define input and output variables
    accepted_in_var = ['pressure', 'traction']
    accepted_out_var = ['displacement']

    @tools.time_initialize
    def __init__(self, parameters):
        super().__init__(parameters)

        if self.version is None:
            raise NotImplementedError(
                'Base class method called, class variable version needs to be set in the derived class')

        # set parameters
        self.settings = parameters['settings']
        self.delta_t = self.settings['delta_t']
        self.timestep_start = self.settings['timestep_start']
        self.save_results = self.settings.get('save_results', 1)  # TODO: implement save_result -> frequency
        self.save_restart = self.settings['save_restart']  # TODO: implement restart
        self.input_file = self.settings['input_file']
        self.cores = self.settings['cores']  # number of CPUs Abaqus has to use
        self.dimensions = self.settings['dimensions']
        self.surfaces = self.settings['surfaces']
        self.port = self.settings.get('port', 40000)
        self.dir_csm = os.path.realpath(self.settings['working_directory'])
        self.path_src = os.path.realpath(os.path.dirname(__file__))
        self.tmp_dir = os.environ.get('TMPDIR', '/tmp')
        self.tmp_dir_unique = os.path.join(self.tmp_dir, f'coconut_{getuser()}_{os.getpid()}_abaqus')

        # initialize coconut messages
        self.coco_messages = tools.CocoMessages(self.dir_csm)
        self.coco_messages.remove_all_messages()

        # initialize variables
        self.env = None
        self.timestep = self.timestep_start
        self.iteration = None
        self.model = None
        self.interface_input = None
        self.interface_output = None
        self.number_of_mps = len(self.surfaces)

        # setting for compilation of solver
        self.compile_clean = self.settings.get('compile_clean', False)
        self.solver = 'AbaqusWrapper'  # name of executable
        self.solver_dir = join(self.path_src, self.solver)
        self.dir_abqw = join(self.dir_csm, self.solver)  # solver Abaqus wrapper directory
        self.compile_log = os.path.join(self.dir_abqw, 'make.log')

        # directories and processes
        self.dir_cse = join(self.dir_csm, 'CSE')  # CSE directory
        if not os.path.exists(self.dir_cse):
            os.mkdir(self.dir_cse)
        if not os.path.exists(self.dir_abqw):
            os.mkdir(self.dir_abqw)
        self.cse_process = None
        self.abaqus_process = None
        self.abaqus_wrapper_process = None

    @tools.time_initialize
    def initialize(self):
        super().initialize()

        # prepare abaqus_v6.env
        with open(join(self.path_src, 'abaqus_v6.env'), 'r') as infile:
            with open(join(self.dir_csm, 'abaqus_v6.env'), 'w') as outfile:
                for line in infile:
                    line = line.replace('|TMP_DIRECTORY_NAME|', self.tmp_dir_unique)
                    if '|' in line:
                        raise ValueError(f'The following line in abaqus_v6.env still contains a \'|\' after '
                                         f'substitution: \n \t{line} \nA parameter was not substituted')
                    outfile.write(line)

        # clean/create tmp directory
        shutil.rmtree(join('/tmp', self.tmp_dir_unique), ignore_errors=True)
        os.mkdir(join('/tmp', self.tmp_dir_unique))

        # verify provided surface names
        self.verify_surface_names()

        # prepare input file for Abaqus
        template_input_file = join(self.dir_csm, self.input_file)
        abaqus_input_file = join(self.dir_csm, 'Abaqus.inp')
        self.prepare_input_file(template_input_file, abaqus_input_file)

        # prepare CSE config xml file
        template_xml_file = join(self.path_src, 'CSE_config.xml')
        cse_config_xml_file = join(self.dir_cse, 'CSE_config.xml')
        self.prepare_xml(template_xml_file, cse_config_xml_file)

        # compile solver AbaqusWrapper if necessary
        self.compile_solver()

        # write input for solver AbaqusWrapper
        solver_input_file = join(self.dir_abqw, 'AbaqusWrapper_input.txt')
        with open(solver_input_file, 'w') as f:
            out = f'dt {self.delta_t}\ndimensions {self.dimensions}\nnumber_of_model_parts ' \
                  f'{len(self.settings["interface_output"])}\nport {self.port}\ndebug {int(self.debug)}\n'
            f.write(out)

        # launch processes
        debug_cmd = f'export ABAQUS_MPF_DIAGNOSTIC_LEVEL=1023 && ' if self.debug else ''

        # launch CSE
        cse_log_file = join(self.dir_cse, 'CSE.log')
        launch_cmd = f'abaqus cse -config {cse_config_xml_file} -listenerPort {self.port} -timeout 86400 &> {cse_log_file}'
        cmd = debug_cmd + launch_cmd
        self.cse_process = subprocess.Popen(cmd, executable='/bin/bash', shell=True, cwd=self.dir_cse, env=self.env)

        # launch Abaqus
        abq_log_file = join(self.dir_csm, 'abaqus.log')
        launch_cmd = f'abaqus job=Abaqus input=Abaqus.inp -cseDirector localhost:{self.port} -timeout 86400 ' \
                     f'interactive  cpus={self.cores} output_precision=full &> {abq_log_file}'
        cmd = debug_cmd + launch_cmd  # cpus={self.cores} output_precision=full interactive
        self.abaqus_process = subprocess.Popen(cmd, executable='/bin/bash', shell=True, cwd=self.dir_csm, env=self.env)

        # TODO: implement check if licensing is ok

        # # check log for completion and or errors
        # subprocess.run(f'tail -n 10 {self.logfile} > Temp_log.coco', shell=True, cwd=self.dir_csm,
        #                executable='/bin/bash', env=self.env)
        # log_tmp = os.path.join(self.dir_csm, 'Temp_log.coco')
        # bool_lic = True
        # with open(log_tmp, 'r') as fp:
        #     for line in fp:
        #         if any(x in line for x in ['Licensing error', 'license error',
        #                                    'Error checking out Abaqus license']):
        #             bool_lic = False
        # if not bool_lic:
        #     tools.print_info('Abaqus licensing error', layout='fail')
        # elif 'COMPLETED' in line:  # check final line for completed
        #     bool_completed = True
        # elif bool_lic:  # final line did not contain 'COMPLETED' but also no licensing error detected
        #     raise RuntimeError(f'Abaqus did not complete, unclassified error, see {self.logfile} for extra '
        #                        f'information')

        # launch AbaqusWrapper
        abqw_log_file = join(self.dir_abqw, f'{self.solver}.log')
        cmd = f'abaqus {join(self.solver_dir, self.solver)} &> {abqw_log_file}'
        self.abaqus_wrapper_process = subprocess.Popen(cmd, executable='/bin/bash', shell=True, cwd=self.dir_abqw,
                                                       env=self.env)

        # pass on process to coco_messages for polling
        self.coco_messages.set_process(self.abaqus_wrapper_process)

        # wait for solver AbaqusWrapper to write mesh data
        self.coco_messages.wait_message('export_mesh_data_ready')

        # create Model
        self.model = data_structure.Model()

        # create input ModelParts (element centroids)
        for mp_id, mp_var in enumerate(self.settings['interface_input']):
            ids_coords = np.loadtxt(join(self.dir_csm, f'initial_element_centroid_coordinates_mp{mp_id}.txt'),
                                    skiprows=2)
            ids_coords = np.hstack((ids_coords, np.zeros((ids_coords.shape[0], 4 - ids_coords.shape[1]))))
            ids, x0, y0, z0 = ids_coords.T
            self.model.create_model_part(mp_var['model_part'], x0, y0, z0, np.round(ids).astype(int))

        # create output ModelParts (nodes)
        for mp_id, mp_var in enumerate(self.settings['interface_output']):
            ids_coords = np.loadtxt(join(self.dir_csm, f'initial_node_coordinates_mp{mp_id}.txt'), skiprows=2)
            ids_coords = np.hstack((ids_coords, np.zeros((ids_coords.shape[0], 4 - ids_coords.shape[1]))))
            ids, x0, y0, z0 = ids_coords.T
            self.model.create_model_part(mp_var['model_part'], x0, y0, z0, np.round(ids).astype(int))

        # create Interfaces
        self.interface_input = data_structure.Interface(self.settings['interface_input'], self.model)
        self.interface_output = data_structure.Interface(self.settings['interface_output'], self.model)

    def initialize_solution_step(self):
        super().initialize_solution_step()

        self.iteration = 0
        self.timestep += 1

        self.coco_messages.send_message('next')
        self.coco_messages.wait_message('next_ready')

    @tools.time_solve_solution_step
    def solve_solution_step(self, interface_input):
        self.iteration += 1

        # store incoming loads
        self.interface_input = interface_input.copy()

        # write loads to a file read by AbaqusWrapper
        comps = {1: '', 2: 'x y', 3: 'x y z'}
        for mp_id, mp_var in enumerate(self.settings['interface_input']):
            mp_name = mp_var['model_part']
            for var in mp_var['variables']:
                file_name = join(self.dir_csm, f'{var}_mp{mp_id}.txt')
                cols = min(data_structure.variables_dimensions[var], self.dimensions)
                data = self.interface_input.get_variable_data(mp_name, var)[:, :cols]
                np.savetxt(file_name, data, header=f'mp{mp_id} {mp_name} {var}: timestep {self.timestep}, iteration '
                                                   f'{self.iteration} #{self.model.get_model_part(mp_name).size}\n'
                                                   f'{comps[cols]}')

        # copy input data for debugging
        if self.debug:
            for mp_id, mp_var in enumerate(self.settings['interface_input']):
                self.copy_for_debugging(join(self.dir_csm, f'pressure_mp{mp_id}.txt'))
                self.copy_for_debugging(join(self.dir_csm, f'traction_mp{mp_id}.txt'))

        # let AbaqusWrapper run, wait for data
        self.coco_messages.send_message('continue')
        self.coco_messages.wait_message('continue_ready')

        # if self.check_coupling_convergence:
        #     # check if Fluent converged after 1 iteration
        #     self.coupling_convergence = self.coco_messages.check_message('solver_converged')
        #     if self.print_coupling_convergence and self.coupling_convergence:
        #         tools.print_info(f'{self.__class__.__name__} converged')

        # read data from AbaqusWrapper
        for mp_id, mp_var in enumerate(self.settings['interface_output']):
            for var in mp_var['variables']:
                file_name = join(self.dir_csm, f'{var}_mp{mp_id}.txt')
                data = np.loadtxt(file_name, skiprows=2)
                data = np.hstack((data, np.zeros((data.shape[0], 3 - data.shape[1]))))
                self.interface_output.set_variable_data(mp_var['model_part'], var, data)

        # copy output data for debugging
        if self.debug:
            for mp_id, mp_var in enumerate(self.settings['interface_output']):
                self.copy_for_debugging(join(self.dir_csm, f'displacement_mp{mp_id}.txt'))

        # return interface_output object
        return self.interface_output

    def finalize_solution_step(self):
        super().finalize_solution_step()

        self.coco_messages.send_message('end_step')
        self.coco_messages.wait_message('end_step_ready')
        # if self.timestep - 1 > self.timestep_start and \
        #         (self.save_restart == 0 or (self.timestep - 1) % self.save_restart != 0):
        #     # no files from previous time step needed for restart
        #     self.remove_files(self.timestep - 1)
        # if self.save_restart < 0 and self.timestep + self.save_restart > self.timestep_start and \
        #         self.timestep % self.save_restart == 0:
        #     # if (self.timestep + self.save_restart < self.timestep_start): don't touch files from previous
        #     # calculation
        #     # files from (self.timestep + self.save_restart) may be removed as new restart files are
        #     # present at current timestep
        #     self.remove_files(self.timestep + self.save_restart)

    @tools.time_save
    def output_solution_step(self):
        super().output_solution_step()
        pass

        # # save if required
        # if (self.save_results != 0 and self.timestep % self.save_results == 0) \
        #         or (self.save_restart != 0 and self.timestep % self.save_restart == 0):
        #     self.coco_messages.send_message('save')
        #     self.coco_messages.wait_message('save_ready')
        #
        # # remove unnecessary files
        # if self.timestep - 1 > self.timestep_start:
        #     self.remove_dat_files(self.timestep - 1)
        #     if self.save_restart < 0 and self.timestep + self.save_restart > self.timestep_start and \
        #             self.timestep % self.save_restart == 0 \
        #             and (self.save_results == 0 or (self.timestep + self.save_restart) % self.save_results != 0):
        #         # new restart file is written (self.timestep % self.save_restart ==0),
        #         # so previous one (at self.timestep + self.save_restart) can be deleted if:
        #         # - save_restart is negative
        #         # - files from a previous calculation are not touched
        #         # - files are not kept for save_results
        #         for extension in ('cas.h5', 'cas', 'dat.h5', 'dat'):
        #             try:
        #                 os.remove(join(self.dir_cfd, f'case_timestep{self.timestep + self.save_restart}.{extension}'))
        #             except OSError:
        #                 continue

    def finalize(self):
        super().finalize()
        self.coco_messages.send_message('stop')
        self.coco_messages.wait_message('stop_ready')
        self.abaqus_wrapper_process.wait()
        self.cse_process.wait()
        self.abaqus_process.wait()

        # remove unnecessary files
        shutil.rmtree(self.tmp_dir_unique, ignore_errors=True)

    def check_software(self):
        """Check whether the software requirements for this wrapper are fulfilled."""
        # Abaqus version
        result = subprocess.run(['abaqus information=release'], shell=True, stdout=subprocess.PIPE, env=self.env)
        if self.version not in str(result.stdout):
            raise RuntimeError(f'Abaqus version {self.version} is required\nCheck if the solver load commands for'
                               f' the "machine_name" are correct in solver_modules.py')

        # compiler
        try:
            result = subprocess.run('g++ --version', shell=True, capture_output=True, env=self.env, check=True, text=True)
            version = re.search(r'g\+\+ \(.*?\) (\d+)\.', result.stdout).group(1)
            if version is None or int(version) < 8:
                warnings.warn(f'The g++ compiler version is\n{result.stdout}\nThis might be insufficient to compile '
                              f'{self.solver}\nCheck {self.compile_log}', category=UserWarning)
        except subprocess.CalledProcessError:
            raise RuntimeError('Intel compiler g++ must be available')

    def verify_surface_names(self):
        subprocess.run(f'abaqus job=datacheck input={self.input_file.removesuffix(".inp")} datacheck', shell=True,
                       cwd=self.dir_csm, env=self.env, check=True)

        data_file = join(self.dir_csm, 'datacheck.dat')
        while not os.path.exists(data_file):  # TODO: implement function to wait for file
            time.sleep(0.5)
        time.sleep(1.0)
        defined_surfaces = set()
        with open(data_file, 'r') as f:
            for line in f:
                defined_surfaces.update(re.findall(r'\s*\*surface, type=ELEMENT, name=(\S+)', line))

        for surface in self.surfaces:
            if surface not in defined_surfaces:
                raise ValueError(f'Surface {surface} is not defined in Abaqus\n'
                                 f'\tAvailable surfaces are:{" ,".join(defined_surfaces)}')

        if not self.debug:  # TODO: fix removal
            to_be_removed_suffix = ['.023', '1.SMABulk', '.cax', '.com', '.log', '.dat', '.mdl', '.msg', '.odb', '.prt',
                                    '.res', '.sim', '.stt']
            for suffix in to_be_removed_suffix:
                try:
                    os.remove(join(self.dir_csm, f'datacheck{suffix}'))
                except FileNotFoundError:
                    pass
                    # print(f'datacheck{suffix} not found')

    def compile_solver(self):
        # compile abaqus wrapper
        source_file = f'{self.solver}.cpp'
        make_cmd = f'abaqus make job={source_file} {"recompile" if self.compile_clean else ""} &> {self.compile_log}'
        try:
            subprocess.run(make_cmd, cwd=self.solver_dir, check=True, shell=True, env=self.env)
        except subprocess.CalledProcessError:
            raise RuntimeError(f'Compilation of {source_file} failed. Check {self.compile_log}')

    def prepare_xml(self, template_xml_file, xml_file):
        tree = ElTr.parse(template_xml_file)
        root = tree.getroot()

        def add_connector(parent_element, name, component_instance, region, in_out_put, variables):
            connector_elem = ElTr.SubElement(parent_element, 'connector', name=name)
            component_instance_elem = ElTr.SubElement(connector_elem, 'componentInstance')
            component_instance_elem.text = component_instance
            variables_elem = ElTr.SubElement(connector_elem, 'variables')
            in_out_put_elem = ElTr.SubElement(variables_elem, in_out_put)
            for var in variables:
                var_elem = ElTr.SubElement(in_out_put_elem, 'variable')
                var_elem.text = f'{var}.{region}'

        def add_connection_set(parent_element, name, connector_1, connector_2):
            connection_set_elem = ElTr.SubElement(parent_element, 'connectionSet', name=name, type='FIELD')
            connection_elem = ElTr.SubElement(connection_set_elem, 'connection', mapOn='NOMAP')
            for connector in (connector_1, connector_2):
                connector_elem = ElTr.SubElement(connection_elem, 'connector')
                connector_elem.text = connector

        def add_to_connection_groups(conn_group_elem, connection_sets):
            for conn_cat_elem in conn_group_elem:
                for connection_set in connection_sets:
                    connection_set_elem = ElTr.SubElement(conn_cat_elem, 'connectionSet')
                    connection_set_elem.text = connection_set

        def replace_parameters_in_xml_attributes(xml_root: ElTr.Element, tag: str, parameter: str, value):
            # replaces all parameters (e.g. |DT|) in elements with provided tag (e.g. constantDt)
            tag_instances = xml_root.findall(f'.//{tag}')
            for tag_instance in tag_instances:
                tag_instance.text = tag_instance.text.replace(parameter, str(value))

        load_vars = ['pressure', 'traction_vector']
        disp_vars = ['displacement']

        connectors_elem = root.find('.//connectors')
        connection_sets_elem = root.find('.//connectionSets')
        connection_groups_elem = root.find('.//connectionGroups')

        # clear elements
        connectors_elem.clear()
        connection_sets_elem.clear()
        for conn_cat_elem in connection_groups_elem:
            while len(conn_cat_elem) > 0:
                conn_cat_elem.remove(conn_cat_elem[0])

        # add elements
        for mp_id in range(self.number_of_mps):
            add_connector(connectors_elem, f'FromAbaqusWrapper{mp_id}', 'AbaqusWrapper', f'wrapperMesh{mp_id}',
                          'output', load_vars)
            add_connector(connectors_elem, f'ToAbaqusWrapper{mp_id}', 'AbaqusWrapper', f'wrapperMesh{mp_id}', 'input',
                          disp_vars)
            add_connector(connectors_elem, f'FromAbaqus{mp_id}', 'Abaqus', self.surfaces[mp_id], 'output', disp_vars)
            add_connector(connectors_elem, f'ToAbaqus{mp_id}', 'Abaqus', self.surfaces[mp_id], 'input', load_vars)
            add_connection_set(connection_sets_elem, f'FromAbaqusWrapperToAbaqus{mp_id}', f'FromAbaqusWrapper{mp_id}',
                               f'ToAbaqus{mp_id}')
            add_connection_set(connection_sets_elem, f'FromAbaqusToAbaqusWrapper{mp_id}', f'FromAbaqus{mp_id}',
                               f'ToAbaqusWrapper{mp_id}')
            add_to_connection_groups(connection_groups_elem,
                                     (f'FromAbaqusWrapperToAbaqus{mp_id}', f'FromAbaqusToAbaqusWrapper{mp_id}'))

        # replace parameters
        replace_parameters_in_xml_attributes(root, 'constantDt', '|DT|', self.delta_t)
        replace_parameters_in_xml_attributes(root, 'duration', '|DURATION|', 'INF')

        # pretty print with indentation
        parsed_xml = minidom.parseString(ElTr.tostring(root, xml_declaration=True, encoding='UTF-8'))
        pretty_xml = parsed_xml.toprettyxml(indent=4 * ' ')
        pretty_xml = '\n'.join(line for line in pretty_xml.split('\n') if line.strip())  # remove empty lines
        with open(xml_file, 'w') as f:
            f.write(pretty_xml)

    def prepare_input_file(self, template_input_file, input_file):
        # write cosimulation settings after step analysis definition
        export_lines = "\n".join([f'{surface}, U' for surface in self.surfaces])
        import_lines = "\n".join([f'{surface}, P, TRVEC' for surface in self.surfaces])
        out = f'**\n' \
              f'** **********************************\n' \
              f'**	     COSIMULATION SETTINGS\n' \
              f'** **********************************\n' \
              f'*co-simulation, name=std_std,program=MULTIPHYSICS\n' \
              f'*co-simulation region, type=SURFACE, export\n' \
              f'{export_lines}\n' \
              f'*co-simulation region, type=SURFACE, import\n' \
              f'{import_lines}\n' \
              f'**\n' \
              f'** **********************************\n'

        # TODO: INITIAL=NO setting?
        with open(template_input_file, 'r') as in_file:
            with open(input_file, 'w') as out_file:
                in_step = False
                analysis_seen = False
                cosimulation_settings_written = False
                incrementation_seen = False
                for line in in_file:
                    if line.startswith('*Step'):
                        in_step = True
                    elif in_step and line[0] == '*' and line[1] != '*' and not analysis_seen:
                        if not (line.lower().startswith('*dynamic') or line.lower().startswith('*static')):
                            warnings.warn(f'Expected *DYNAMIC or *STATIC instead of\n\t{line}\nCheck resulting input '
                                          f'file {input_file} to see if insertions were done correctly',
                                          category=UserWarning)
                        analysis_seen = True
                    if in_step and analysis_seen and not line.startswith('*') and not incrementation_seen:
                        # on data line for time increment
                        incr_info = [inc.strip() for inc in line.split(',')]
                        if not np.isclose(float(incr_info[0]), self.delta_t, atol=0):
                            warnings.warn(f'Suggested time increment {float(incr_info[0])} in input file not used '
                                          f'because increment is set by CSE: delta_t = {self.delta_t} ',
                                          category=UserWarning)
                            incr_info[0] = str(self.delta_t)
                        incr_info[1] = '999999'
                        line_new = ', '.join(incr_info) + '\n'
                        out_file.write(line_new)
                        if not cosimulation_settings_written:
                            out_file.write(out)
                            cosimulation_settings_written = True
                        incrementation_seen = True
                    else:
                        out_file.write(line)

    def copy_for_debugging(self, path):
        new_path = f'{path[:-len(".txt")]}_ts{self.timestep}_it{self.iteration}.txt'
        shutil.copy2(path, new_path)
