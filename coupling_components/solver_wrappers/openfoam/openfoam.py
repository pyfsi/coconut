from coconut import data_structure
from coconut.coupling_components.component import Component
from coconut.data_structure.interface import Interface
from coconut import tools
from coconut.coupling_components.solver_wrappers.openfoam import openfoam_io as of_io
import matplotlib.pyplot as plt

import copy
import numpy as np
import os
import shutil
import time
import subprocess
import re
from glob import glob


def create(parameters):
    return SolverWrapperOpenFOAM(parameters)


class SolverWrapperOpenFOAM(Component):
    version = None  # OpenFOAM version with dot, e.g. 4.1 , set in sub-class

    @tools.time_initialize
    def __init__(self, parameters):
        super().__init__()

        if self.version is None:
            raise NotImplementedError(
                'Base class method called, class variable version needs to be set in the derived class')

        # settings
        self.settings = parameters['settings']
        self.working_directory = self.settings['working_directory']
        self.env = None  # environment in which correct version of OpenFOAM software is available, set in sub-class

        # adapted application from openfoam ('coconut_<application name>')
        self.application = self.settings['application']
        self.delta_t = self.settings['delta_t']
        self.number_of_timesteps = self.settings['number_of_timesteps']
        self.time_precision = self.settings['time_precision']
        self.start_time = self.settings['timestep_start'] * self.delta_t
        self.timestep = self.physical_time = self.iteration = self.prev_timestamp = self.cur_timestamp = None
        self.save_restart = self.settings['save_restart']
        self.openfoam_process = None
        self.write_interval = self.write_precision = None
        self.adjustFSITimeStep = None
        self.velocityList=[]

        # boundary_names is the set of boundaries in OpenFoam used for coupling
        self.boundary_names = self.settings['boundary_names']
        self.moving_rigid_bodies_names = self.settings['moving_rigid_bodies_names']
        self.cores = None
        self.model = None
        self.interface_input = None
        self.interface_output = None

        # set on True to save copy of input and output files in every iteration
        self.debug = self.settings.get('debug', False)

        # set on True if you want to clean the adapted application and compile.
        self.compile_clean = self.settings.get('compile_clean', False)

        # remove possible CoCoNuT-message from previous interrupt
        self.remove_all_messages()

        # variables for kinematic pressure and traction conversion
        self.density_for_pressure = None
        self.density_for_traction = None
        self.wall_shear_stress_variable = None
        self.wall_shear_stress_function_object_library = None

        # time
        self.init_time = self.init_time
        self.run_time = 0.0

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
                    f'In the parameters:\n{self.settings}\n key "parallel" is set to {True} but {file_name} '
                    f'does not exist')
            else:
                with open(file_name, 'r') as file:
                    decomposedict_string = file.read()
                self.cores = of_io.get_int(input_string=decomposedict_string, keyword='numberOfSubdomains')

        # based on solver application, set conversion settings from kinematic to static pressure/shear stress
        # typically: incompressible solver, pressure and wallShearStress are kinematic -> multiply with fluid density
        #            compressible solver, pressure and wallShearStress are not kinematic -> do nothing
        kinematic_conversion_dict = {
            'coconut_cavitatingFoam': {
                'wall_shear_stress_variable': 'rhoWallShearStress'
            },
            'coconut_mod_cavitatingFoam': {
                'wall_shear_stress_variable': 'rhoWallShearStress'
            },
            'coconut_mod_cavitatingFoam_therm': {
                'wall_shear_stress_variable': 'rhoWallShearStress'
            },
            'coconut_interFoam': {
                'wall_shear_stress_variable': 'rhoWallShearStress'
            },
            'coconut_pimpleFoam': {
                'density_for_pressure': 'look up',
                'density_for_traction': 'look up',
            }
        }
        if self.application not in kinematic_conversion_dict:
            available_applications = ''
            for key in kinematic_conversion_dict:
                available_applications += f'\n\t{key}'
            raise ValueError(f'{self.application} is not included in the kinematic_conversion_dict '
                             f'used for treatment of (kinematic) pressure and traction\n'
                             f'Add the solver to this dictionary '
                             f'or use one of the existing solvers:{available_applications}')
        else:
            kinematic_conversion = kinematic_conversion_dict[self.application]
            self.density_for_pressure = 1.0 if 'density_for_pressure' not in kinematic_conversion \
                else self.settings['density']  # default density is 1
            self.density_for_traction = 1.0 if 'density_for_traction' not in kinematic_conversion \
                else self.settings['density']  # default density is 1
            self.wall_shear_stress_variable = kinematic_conversion.get(
                'wall_shear_stress_variable', 'wallShearStress')  # default shear stress variable is wallShearStress

        # modify controlDict file to add pressure and wall shear stress functionObjects for all the boundaries in
        # self.settings["boundary_names"]
        self.read_modify_controldict()

        if self.save_restart % self.write_interval:
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
            self.update_axial_velocity = np.ones(2*nfaces +2)

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

        if self.settings['timeVaryingMappedFixedValue']:
            for boundary in self.boundary_names:
                mp_name = f'{boundary}_input'
                mp = self.model.get_model_part(mp_name)
                x0, y0, z0 = mp.x0, mp.y0, mp.z0

                x = np.zeros(x0.size)
                y = np.zeros(x0.size)
                z = np.zeros(x0.size)
                index = int(len(x)/2)

                j = 0
                for i in range(len(x)):
                    if z0[i] < 0:
                        x[j] = x0[i]
                        y[j] = y0[i]
                        z[j] = z0[i]

                    else:
                        x[j + index] = x0[i]
                        y[j + index] = y0[i]
                        z[j + index] = z0[i]
                        j += 1

                boundary_data_path = os.path.join(self.working_directory, 'constant/boundaryData')
                if not os.path.exists(boundary_data_path):
                    os.mkdir(boundary_data_path)
                boundary_path = os.path.join(boundary_data_path, boundary)
                shutil.rmtree(boundary_path, ignore_errors=True)
                os.mkdir(boundary_path)
                data_folder = os.path.join(boundary_path, '0')
                os.mkdir(data_folder)

                with open(os.path.join(boundary_path, 'points'), 'w') as f:
                    f.write('(\n')
                    for point in range(x.size):
                        f.write(f'({x[point]} {y[point]} {z[point]})\n')
                    f.write(')')

                with open(os.path.join(data_folder, 'U'), 'w') as h:
                    h.write(f'{x.size}\n')
                    h.write('(\n')
                    for i in range(x.size):
                        h.write(f'({1} {0} {0} )\n')
                    h.write(')')

        # define timestep and physical time
        self.timestep = 0
        self.physical_time = self.start_time

        # copy zero folder to folder with correctly named timeformat
        if self.start_time == 0:
            timestamp = '{:.{}f}'.format(self.physical_time, self.time_precision)
            path_orig = os.path.join(self.working_directory, '0')
            path_new = os.path.join(self.working_directory, timestamp)
            shutil.rmtree(path_new, ignore_errors=True)
            shutil.copytree(path_orig, path_new)
            if self.settings['moving_rigid_body']:
                for rigid in self.moving_rigid_bodies_names:
                    self.radial_displacement = self.settings['radial_displacement_rigid_body']
                    self.number_of_timeIncrements = self.settings['number_of_time_increments_rigid_body']  # number of timeIncrements is the value how many times a displacement of the rigid body is done.
                    self.radius = self.settings['radius_wire'],
                    self.start_increment = self.settings['start_increment']
                    # rigid_body_filename = os.path.join( self.working_directory,'constant/polyMesh/boundary')
                    for rigid in self.moving_rigid_bodies_names:
                        # with open (rigid_body_filename,'r') as rigidBody_file:
                        # rigidBody_string = rigidBody_file.read()
                        # rigidBody_dict = of_io.get_dict(input_string=rigidBody_string, keyword=rigid)
                        # get point ids and coordinates for all the faces in the boundary
                        node_ids, node_coords = of_io.get_boundary_points(case_directory=self.working_directory,
                                                                          time_folder='0',
                                                                          boundary_name=rigid)

                        rigid_model = data_structure.Model()

                        # create input model part
                        rigid_model.create_model_part(f'{rigid}', node_coords[:, 0], node_coords[:, 1],
                                                      node_coords[:, 2],
                                                      node_ids)

                        mp_name = f'{rigid}'
                        mp = rigid_model.get_model_part(mp_name)
                        x0, y0, z0 = mp.x0, mp.y0, mp.z0

                        x0[0], x0[1] = x0[1], x0[0]
                        y0[0], y0[1] = y0[1], y0[0]
                        z0[0], z0[1] = z0[1], z0[0]

                        self.x_rigid = np.zeros(x0.size)
                        self.y_rigid = np.zeros(x0.size)
                        self.z_rigid = np.zeros(x0.size)
                        index = int(len(self.x_rigid) / 2)

                        j = 0
                        for i in range(len(self.x_rigid)):
                            if z0[i] > 0:
                                self.x_rigid[j + index] = x0[i]
                                self.y_rigid[j + index] = y0[i]
                                self.z_rigid[j + index] = z0[i]

                            else:
                                self.x_rigid[j] = x0[i]
                                self.y_rigid[j] = y0[i]
                                self.z_rigid[j] = z0[i]
                                j += 1

                        rigid_data_path = os.path.join(self.working_directory, 'constant/boundaryData')
                        if not os.path.exists(rigid_data_path):
                            os.mkdir(rigid_data_path)
                        rigid_path = os.path.join(rigid_data_path, rigid)
                        shutil.rmtree(rigid_path, ignore_errors=True)
                        os.mkdir(rigid_path)
                        data_folder = os.path.join(rigid_path, '0')
                        os.mkdir(data_folder)

                        with open(os.path.join(rigid_path, 'points'), 'w') as f:
                            f.write('(\n')
                            for point in range(self.x_rigid.size):
                                f.write(f'({self.x_rigid[point]} {self.y_rigid[point]} {self.z_rigid[point]})\n')
                            f.write(')')

                        with open(os.path.join(data_folder, 'pointDisplacement'), 'w') as h:
                            h.write(f'{self.x_rigid.size}\n')
                            h.write('(\n')
                            for i in range(self.x_rigid.size):
                                h.write(f'({0} {0} {0} )\n')
                            h.write(')')
                        timestamp = '{:.{}f}'.format(self.physical_time, self.time_precision)
                        path_orig_rigidData = os.path.join(self.working_directory, 'constant/boundaryData', rigid, '0')
                        path_new_rigidData = os.path.join(self.working_directory,'constant/boundaryData', rigid, timestamp )
                        shutil.rmtree(path_new_rigidData, ignore_errors=True)
                        shutil.copytree(path_orig_rigidData, path_new_rigidData)

                       #Make a cosine function for softer displacement increments during die diameter decrease. The cosine is added after the
                       # linear decreas of the die diameter. Total displacement is represented by the parameter
                        rico = (self.radius[0] - (self.radius[0]-self.radial_displacement))/ (self.delta_t * self.number_of_timeIncrements)
                        x = np.linspace(0,self.delta_t*self.number_of_timeIncrements,self.number_of_timeIncrements)
                        y = -1*rico*x
                        displacement_increment = np.zeros([len(x) + self.number_of_timeIncrements//10])
                        q = self.start_increment*self.delta_t + np.linspace(0, self.delta_t * self.number_of_timeIncrements + (self.number_of_timeIncrements//10) * self.delta_t, len(displacement_increment))
                        p = self.start_increment*self.delta_t + x

                        delta_displacement = self.radial_displacement / self.number_of_timeIncrements
                        total_delta = self.radial_displacement
                        for g in range(len(displacement_increment)):
                            if g < self.number_of_timeIncrements:
                                displacement_increment[g] = y[g]
                            else:
                                new_delta = delta_displacement/1.09
                                delta_displacement = new_delta
                                total_delta +=new_delta
                                displacement_increment[g] = -total_delta

                        for i in range(self.number_of_timeIncrements + self.number_of_timeIncrements//10):
                            # print("i")
                            # print(i)
                            if self.adjustFSITimeStep:
                            #add variable j to achieve the correct time for adjustable timesteps in order to decrease the die gap on the correct timestep
                                j = 0
                                while j < self.start_increment:
                                    j+=1
                                    if j < 51:
                                        t = self.delta_t * (j + 1)
                                    else:
                                        delta_time = self.delta_t * 0.08 ** (j - 50)
                                        if delta_time < 1e-8:
                                            delta_time = 1e-8
                                        t += delta_time
                                time = t + delta_time * i
                            else:
                                time = self.delta_t * (i + 1 + self.start_increment)
                            # print("time")
                            # print(time)

                            timestamp = '{:.{}f}'.format(time, self.time_precision)
                            new_path_rigidData = os.path.join(self.working_directory, 'constant/boundaryData', rigid, timestamp)
                            shutil.rmtree(new_path_rigidData, ignore_errors=True)
                            shutil.copytree(path_new_rigidData, new_path_rigidData)

                            deltaY = np.zeros(len(self.y_rigid))
                            deltaZ = np.zeros(len(self.z_rigid))

                            for j in range(len(self.z_rigid)):
                               deltaY[j] = displacement_increment[i] * np.cos(2.5 * np.pi / 180)
                               if self.z_rigid[j] < 0:
                                   deltaZ[j] = -displacement_increment[i] * np.sin(2.5 * np.pi / 180)
                               else:
                                   deltaZ[j] = displacement_increment[i] * np.sin(2.5 * np.pi / 180)

                            with open(os.path.join(new_path_rigidData,'pointDisplacement'), 'w') as h:
                                h.write(f'{self.x_rigid.size}\n')
                                h.write('(\n')
                                for k in range(self.x_rigid.size):
                                    h.write(f'({0} {deltaY[k]} {deltaZ[k]})\n')
                                h.write(')')

            # plt.plot(q, displacement_increment)
            # plt.show()

            if self.settings['timeVaryingMappedFixedValue']:
                for boundary in self.boundary_names:
                    timestamp = '{:.{}f}'.format(self.physical_time, self.time_precision)
                    path_orig_boundaryData = os.path.join(self.working_directory, 'constant/boundaryData', boundary, '0')
                    path_new_boundaryData = os.path.join(self.working_directory,'constant/boundaryData', boundary, timestamp )
                    shutil.rmtree(path_new_boundaryData, ignore_errors=True)
                    shutil.copytree(path_orig_boundaryData, path_new_boundaryData)
                    for i in range(self.number_of_timesteps + 1):
                        if self.adjustFSITimeStep:
                            if i < 51:
                                time = self.delta_t * (i + 1)
                            else:
                                delta_time = self.delta_t*0.08**(i-50)
                                if delta_time < 1e-8:
                                    delta_time = 1e-8
                                time += delta_time
                        else:
                            time = self.delta_t * (i + 1)
                        timestamp = '{:.{}f}'.format(time, self.time_precision)
                        new_path_boundaryData = os.path.join(self.working_directory, 'constant/boundaryData', boundary,
                                                                 timestamp)
                        shutil.copytree(path_orig_boundaryData, new_path_boundaryData)


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

                    if self.settings['timeVaryingMappedFixedValue']:

                        path_working_directory_processor = os.path.join(self.working_directory, f'processor{p}')

                        self.write_cell_centres_parallel_timeVaryingMappedCoupledVelocity(path_working_directory_processor)

                        boundary_filename = os.path.join(self.working_directory, f'processor{p}/constant/polyMesh/boundary')
                        for boundary in self.boundary_names:
                            with open(boundary_filename, 'r') as boundary_file:
                                boundary_file_string = boundary_file.read()
                            boundary_dict = of_io.get_dict(input_string=boundary_file_string, keyword=boundary)
                            # get point ids and coordinates for all the faces in the boundary
                            # node_ids, node_coords = of_io.get_boundary_points(case_directory=self.working_directory/ f'processor{p}',
                            #                                                   time_folder='0',
                            #                                                   boundary_name=boundary)
                            nfaces = of_io.get_int(input_string=boundary_dict, keyword='nFaces')

                        x0, y0, z0 = self.read_face_centres_parallel_timeVaryingMappedCoupledVelocity(boundary,nfaces,p)

                        x_p = np.zeros(2* x0.size)
                        y_p = np.zeros(2* x0.size)
                        z_p = np.zeros(2* x0.size)
                        index = int(len(x_p) / 2)


                        j = 0
                        for i in range(len(x_p)):
                            if i < len(x0):
                                x_p[j] = x0[i]
                                y_p[j] = y0[i] * np.cos(np.pi/180)
                                z_p[j] = - y0[i] * np.sin(np.pi/180)

                            else:
                                x_p[j] = x0[i-index]
                                y_p[j] = y0[i-index] * np.cos(np.pi/180)
                                z_p[j] = y0[i-index] * np.sin(np.pi/180)
                            j += 1

                        boundary_data_path = os.path.join(self.working_directory, f'processor{p}/constant/boundaryData')
                        if not os.path.exists(boundary_data_path):
                            os.mkdir(boundary_data_path)
                        boundary_path = os.path.join(boundary_data_path, boundary)
                        shutil.rmtree(boundary_path, ignore_errors=True)
                        os.mkdir(boundary_path)
                        data_folder = os.path.join(boundary_path, '0')
                        os.mkdir(data_folder)

                        with open(os.path.join(boundary_path, 'points'), 'w') as f:
                            f.write('(\n')
                            for point in range(x_p.size):
                                f.write(f'({x_p[point]} {y_p[point]} {z_p[point]})\n')
                            f.write(')')

                        with open(os.path.join(data_folder, 'U'), 'w') as h:
                            h.write(f'{x_p.size}\n')
                            h.write('(\n')
                            for i in range(x_p.size):
                                h.write(f'({1} {0} {0} )\n')
                            h.write(')')

                        timestamp = '{:.{}f}'.format(self.physical_time, self.time_precision)
                        path_orig_boundaryData = os.path.join(self.working_directory,
                                                              f'processor{p}/constant/boundaryData', boundary, '0')
                        # os.makedirs(path_orig_boundaryData)
                        path_new_boundaryData = os.path.join(self.working_directory,
                                                             f'processor{p}/constant/boundaryData', boundary,
                                                             timestamp)
                        shutil.rmtree(path_new_boundaryData, ignore_errors=True)
                        shutil.copytree(path_orig_boundaryData, path_new_boundaryData)
                        for i in range(self.number_of_timesteps + 1):
                            if self.adjustFSITimeStep:
                                if i < 51:
                                    time = self.delta_t * (i + 1)
                                else:
                                    delta_time = self.delta_t * 0.08 ** (i - 50)
                                    if delta_time < 1e-8:
                                        delta_time = 1e-8
                                    time += delta_time
                            else:
                                time = self.delta_t * (i + 1)
                            timestamp = '{:.{}f}'.format(time, self.time_precision)
                            new_path_boundaryData = os.path.join(self.working_directory,
                                                                 f'processor{p}/constant/boundaryData',
                                                                 boundary,
                                                                 timestamp)
                            shutil.copytree(path_orig_boundaryData, new_path_boundaryData)

        # starting the OpenFOAM infinite loop for coupling!
        if not self.settings['parallel']:
            cmd = self.application + '&> log.' + self.application
        else:
            cmd = 'mpirun -np ' + str(self.cores) + ' ' + self.application + ' -parallel &> log.' + self.application

        self.openfoam_process = subprocess.Popen(cmd, cwd=self.working_directory, shell=True, env=self.env)

    def initialize_solution_step(self):
        super().initialize_solution_step()

        timestamp = '{:.{}f}'.format(self.physical_time, self.time_precision)

        if self.settings['timeVaryingMappedFixedValue']:
            for boundary in self.boundary_names:
                path_boundaryData = os.path.join(self.working_directory, 'constant/boundaryData', boundary, timestamp)
                if self.cores > 1 or self.physical_time == 0:
                    os.makedirs(path_boundaryData, exist_ok=True)

        # prepare new time step folder and reset the number of iterations
        self.timestep += 1
        self.iteration = 0

        if self.adjustFSITimeStep:
            if self.timestep < 52:
                self.physical_time += self.delta_t
            else:
                delta_time = self.delta_t * 0.08 ** (self.timestep - 51)
                if delta_time < 1e-8:
                    delta_time = 1e-8
                self.physical_time += delta_time
        else:
            self.physical_time += self.delta_t

        self.prev_timestamp = timestamp
        self.cur_timestamp = f'{self.physical_time:.{self.time_precision}f}'

        if not self.settings['parallel']:  # if serial
            new_path = os.path.join(self.working_directory, self.cur_timestamp)
            if os.path.isdir(new_path):
                tools.print_info(f'Overwrite existing time step folder: {new_path}', layout='warning')
                subprocess.check_call(f'rm -rf {new_path}', shell=True)

            if self.settings['timeVaryingMappedFixedValue']:
                for boundary in self.boundary_names:
                    new_path_boundaryData = os.path.join(self.working_directory, 'constant/boundaryData', boundary,
                                                         self.cur_timestamp)
                if os.path.isdir(new_path):
                    tools.print_info(f'Overwrite existing time step folder: {new_path}', layout='warning')
                    subprocess.check_call(f'rm -rf {new_path}', shell=True)
                if os.path.isdir(new_path_boundaryData):
                    tools.print_info(f'Overwrite existing time step folder: {new_path_boundaryData}', layout='warning')
                    subprocess.check_call(f'rm -rf {new_path_boundaryData}', shell=True)
                subprocess.check_call(f'mkdir -p {new_path}', shell=True)
                subprocess.check_call(f'mkdir -p {new_path_boundaryData}', shell=True)
        else:
            for i in np.arange(self.cores):
                new_path = os.path.join(self.working_directory, 'processor' + str(i), self.cur_timestamp)
                if os.path.isdir(new_path):
                    if i == 0:
                        tools.print_info(f'Overwrite existing time step folder: {new_path}', layout='warning')
                    subprocess.check_call(f'rm -rf {new_path}', shell=True)

                if self.settings['timeVaryingMappedFixedValue']:
                    for boundary in self.boundary_names:
                        new_path_boundaryData = os.path.join(self.working_directory, f'processor{i}/constant/boundaryData', boundary,
                                                             self.cur_timestamp)
                    if os.path.isdir(new_path):
                        if i == 0:
                            tools.print_info(f'Overwrite existing time step folder: {new_path}', layout='warning')
                        subprocess.check_call(f'rm -rf {new_path}', shell=True)
                    if os.path.isdir(new_path_boundaryData):
                        if i == 0:
                            tools.print_info(f'Overwrite existing time step folder: {new_path_boundaryData}',
                                         layout='warning')
                        subprocess.check_call(f'rm -rf {new_path_boundaryData}', shell=True)
                    subprocess.check_call(f'mkdir -p {new_path}', shell=True)
                    subprocess.check_call(f'mkdir -p {new_path_boundaryData}', shell=True)


        self.send_message('next')
        self.wait_message('next_ready')

    @tools.time_solve_solution_step
    def solve_solution_step(self, interface_input):
        self.iteration += 1

        # store incoming displacements
        self.interface_input.set_interface_data(interface_input.get_interface_data())

        # write interface data to OpenFOAM-file
        self.write_node_input()

        # copy input data for debugging
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

        self.send_message('continue')
        self.wait_message('continue_ready')

        # read data from OpenFOAM
        self.read_node_output()

        # copy output data for debugging
        if self.debug:
            for boundary in self.boundary_names:
                post_process_time_folder = os.path.join(self.working_directory,
                                                        f'postProcessing/coconut_{boundary}/surface',
                                                        self.cur_timestamp)
                for file in (f'{self.wall_shear_stress_variable}_patch_{boundary}.raw', f'p_patch_{boundary}.raw'):
                    path = os.path.join(post_process_time_folder, file)
                    iter_path = os.path.join(post_process_time_folder, f'{file[:-4]}_{self.iteration}.raw')
                    shutil.copy(path, iter_path)

        # return interface_output object
        return self.interface_output

    def finalize_solution_step(self):
        super().finalize_solution_step()

        # if os.path.exists(self.file_path):
        #     print("REMOVE")
        #     os.remove(self.file_path)
        #if debugging is not activated, to save memory the post process files are removed after each time step and the time folder and boundaryData folders are saved after each 10 time steps.
        # This has no influence on the solution of the calculation
        if not self.debug:
            for boundary in self.boundary_names:
                post_process_time_folder = os.path.join(self.working_directory,
                                                        f'postProcessing/coconut_{boundary}/surface',
                                                        self.prev_timestamp)
                post_process_time_folder_pressure = os.path.join(self.working_directory,
                                                                 f'postProcessing/PRESSURE_{boundary}/surface',
                                                                 self.prev_timestamp)
                post_process_time_folder_traction = os.path.join(self.working_directory,
                                                                 f'postProcessing/TRACTION_{boundary}/surface',
                                                                 self.prev_timestamp)

                shutil.rmtree(post_process_time_folder)
                shutil.rmtree(post_process_time_folder_pressure)
                shutil.rmtree(post_process_time_folder_traction)



                if (self.timestep % 10 != 1):
                    if self.settings['timeVaryingMappedFixedValue']:
                        prev_directory_boundaryData_coupledVelocity_folder = os.path.join(self.working_directory,
                                                                                          f'constant/boundaryData/{boundary}',
                                                                                          self.prev_timestamp)
                        shutil.rmtree(prev_directory_boundaryData_coupledVelocity_folder)
                    if self.settings[
                        'moving_rigid_body'] and self.start_increment < self.timestep - 1 and self.timestep - 1 < self.start_increment + self.number_of_timeIncrements + self.number_of_timeIncrements // 10:
                        for body in self.moving_rigid_bodies_names:
                            prev_directory_boundaryData_movingBody_folder = os.path.join(self.working_directory,
                                                                                         f'constant/boundaryData/{body}',
                                                                                         self.prev_timestamp)
                            shutil.rmtree(prev_directory_boundaryData_movingBody_folder)

                    if not self.settings['parallel']:
                        prev_directory_folder = os.path.join(self.working_directory, self.prev_timestamp)
                        shutil.rmtree(prev_directory_folder)

                    else:
                        for proc in range(self.cores):
                           prev_directory_folder_proc = os.path.join(self.working_directory, f'processor{proc}',self.prev_timestamp)
                           shutil.rmtree(prev_directory_folder_proc)

                           if self.settings['timeVaryingMappedFixedValue']:
                               prev_directory_boundaryData_coupledVelocity_folder_proc = os.path.join(self.working_directory,f'processor{proc}/constant/boundaryData/{boundary}',self.prev_timestamp)
                               shutil.rmtree(prev_directory_boundaryData_coupledVelocity_folder_proc)

        if not (self.timestep % self.write_interval):
            self.send_message('save')
            self.wait_message('save_ready')

        if self.residual_variables is not None:
            self.write_of_residuals()

    def finalize(self):
        super().finalize()

        self.send_message('stop')
        self.wait_message('stop_ready')
        self.openfoam_process.wait()
        if not self.debug:
            files_to_delete = glob(os.path.join(self.working_directory, 'constant/pointDisplacementTmp*')) + glob(
                os.path.join(self.working_directory, 'processor*/constant/pointDisplacementTmp*'))
            for filepath in files_to_delete:
                os.remove(filepath)

            for boundary in self.boundary_names:
                post_process_folder = os.path.join(self.working_directory, f'postProcessing/coconut_{boundary}')
                shutil.rmtree(post_process_folder, ignore_errors=True)

    def get_interface_input(self):
        return self.interface_input

    def get_interface_output(self):
        return self.interface_output

    def compile_adapted_openfoam_solver(self):
        # compile openfoam adapted solver
        solver_dir = os.path.join(os.path.dirname(__file__), f'v{self.version.replace(".", "")}', self.application)
        boundary_cond_dir = os.path.join(os.path.dirname(__file__), f'v{self.version.replace(".", "")}', 'coconut_src',
                                         'boundaryConditions')
        try:
            if self.compile_clean:
                subprocess.check_call(f'wclean {solver_dir} && wmake {solver_dir} &> log.wmake',
                                      cwd=self.working_directory, shell=True,
                                      env=self.env)
                subprocess.check_call(
                    f'wclean {boundary_cond_dir} && wmake libso {boundary_cond_dir} &> log.wmake_libso',
                    cwd=self.working_directory, shell=True,
                    env=self.env)
            else:
                subprocess.check_call(f'wmake {solver_dir} &> log.wmake', cwd=self.working_directory, shell=True,
                                      env=self.env)
                subprocess.check_call(f'wmake libso {boundary_cond_dir} &> log.wmake_libso',
                                      cwd=self.working_directory, shell=True,
                                      env=self.env)
        except subprocess.CalledProcessError:
            raise RuntimeError(
                f'Compilation of {self.application} or coconut_src failed. Check {os.path.join(self.working_directory, "log.wmake or CFD/log.wmake_libso")}')

    def write_cell_centres(self):
        raise NotImplementedError('Base class method is called, should be implemented in derived class')

    def write_cell_centres_parallel_timeVaryingMappedCoupledVelocity(self):
        raise NotImplementedError('Base class method is called, should be implemented in derived class')

    def read_face_centres(self, boundary_name, nfaces):
        raise NotImplementedError('Base class method is called, should be implemented in derived class')

    def read_face_centres_parallel_timeVaryingMappedCoupledVelocity(self, boundary_name, nfaces):
        raise NotImplementedError('Base class method is called, should be implemented in derived class')

    def delete_prev_iter_output(self):
        # pressure and wall shear stress files are removed to avoid openfoam to append data in the new iteration
        for boundary in self.boundary_names:
            # specify location of pressure and traction
            post_process_time_folder = os.path.join(self.working_directory,
                                                    f'postProcessing/coconut_{boundary}/surface', self.cur_timestamp)
            for file in (f'{self.wall_shear_stress_variable}_patch_{boundary}.raw', f'p_patch_{boundary}.raw'):
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
            x0, y0, z0 = mp.x0, mp.y0, mp.z0
            nfaces = mp.size
            post_process_time_folder = os.path.join(self.working_directory,
                                                    f'postProcessing/coconut_{boundary}/surface', self.cur_timestamp)
            wss_filename = os.path.join(post_process_time_folder, f'{self.wall_shear_stress_variable}_patch_{boundary}.raw')
            pres_filepath = os.path.join(post_process_time_folder, f'p_patch_{boundary}.raw')

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

            wall_shear_stress[pos_list,] = wss_tmp[:, ]
            pressure[pos_list, 0] = pres_tmp

            self.interface_output.set_variable_data(mp_name, 'traction', -wall_shear_stress * self.density_for_traction)
            self.interface_output.set_variable_data(mp_name, 'pressure', pressure * self.density_for_pressure)
            #plt.plot(x0,pressure)
            #plt.show()

    # noinspection PyMethodMayBeStatic
    def write_footer(self, file_name):
        # write OpenFOAM-footer at the end of file
        with open(file_name, 'a') as f:
            f.write('\n// ************************************************************************* //\n')

    def write_node_input(self):
        """
        creates pointDisplacementNext for supplying the displacement field in the FSI coupling. This file is created
        from 0/pointDisplacement. The boundary field for boundaries participating in the FSI coupling is modified to
        supply the boundary displacement field from structural solver. If the OpenFOAM solver is run in parallel,
        the field is subsequently decomposed using the command: decomposePar.
       :return:
       """
        for boundary in self.boundary_names:
            mp_name = f'{boundary}_input'
            displacement = self.interface_input.get_variable_data(mp_name, 'displacement')

            if self.settings['timeVaryingMappedFixedValue']:
                # cwd = os.getcwd()
                # CSM_directory = cwd + '/CSM'
                # self.file_path = CSM_directory + '/update_velocity'
                if self.timestep <= 10:
                    displacement[:, 0] = 1
                    velocity = displacement[:, 0]
                    # speed =velocity
                else:
                    velocity = displacement[:, 0] * 1e5
                velocity = copy.deepcopy(velocity)
                # self.velocityList.append(velocity)
                # if self.timestep >1:
                #     if os.path.exists(self.file_path):
                #         speed = self.velocityList[-1]
                #         # print("Update_velocity")
                #     else:
                #         speed = self.velocityList[-2]
                #         self.velocityList[-1] =speed
                #         # print("No_update_velocity")
                #         self.velocityList.pop(-2)

                # speed = copy.deepcopy(speed)

                data_folder_home = os.path.join(self.working_directory, 'constant/boundaryData', boundary,
                                                self.cur_timestamp)
                velocity_input = np.zeros((len(velocity), 3))
                index = int(len(velocity) / 2)

                j = 0
                for i in range(len(velocity)):
                    if i % 2 == 0:
                        velocity_input[j, 0] = velocity[i]
                    else:
                        velocity_input[j + index, 0] = velocity[i]
                        j += 1

                with open(os.path.join(data_folder_home, 'U'), 'w') as f:
                    f.write(f'{velocity_input.shape[0]}\n')
                    f.write('(\n')
                    for i in range(velocity_input.shape[0]):
                        f.write(f'({velocity_input[i, 0]} {velocity_input[i, 1]} {velocity_input[i, 2]})\n')
                    f.write(')')

                    # for timeLabel in range(self.timestep, self.number_of_timesteps ):
                    #     time = (self.timeList[timeLabel])
                    #     timestamp = '{:.{}f}'.format(time, self.time_precision)
                    #     print("timeStamp")
                    #     print(timestamp)
                    #     updated_boundaryData = os.path.join(self.working_directory, 'constant/boundaryData',
                    #                                          boundary,
                    #                                          timestamp)
                    #     shutil.copytree(data_folder_home, updated_boundaryData)

            displacement[:, 0] = 0

        if not self.settings['parallel']:
            pointdisp_filename_ref = os.path.join(self.working_directory, '0/pointDisplacement')
            pointdisp_filename = os.path.join(self.working_directory, 'constant/pointDisplacementTmp')

            with open(pointdisp_filename_ref, 'r') as ref_file:
                pointdisp_string = ref_file.read()

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

                seq = self.mp_in_decompose_seq_dict[mp_name][proc]

                data_folder = os.path.join(self.working_directory, f'processor{proc}', 'constant/boundaryData',
                                           boundary, self.cur_timestamp)

                if self.settings['timeVaryingMappedFixedValue']:

                    velocity_input = np.zeros((len(velocity[seq]), 3))
                    index = int(len(velocity[seq]) / 2)

                    j = 0
                    for i in range(len(velocity[seq])):
                        if i % 2 == 0:
                            velocity_input[j, 0] = velocity[seq][i]
                        else:
                            velocity_input[j + index, 0] = velocity[seq][i]
                            j += 1

                    with open(os.path.join(data_folder, 'U'), 'w') as f:
                        f.write(f'{velocity_input.shape[0]}\n')
                        f.write('(\n')
                        for i in range(velocity_input.shape[0]):
                            f.write(f'({velocity_input[i, 0]} {velocity_input[i, 1]} {velocity_input[i, 2]})\n')
                        f.write(')')

                if seq is not None:
                    boundary_dict = of_io.get_dict(input_string=pointdisp_string, keyword=boundary)
                    boundary_dict_new = of_io.update_vector_array_dict(dict_string=boundary_dict,
                                                                       vector_array=displacement[seq])
                    pointdisp_string = pointdisp_string.replace(boundary_dict, boundary_dict_new)

                with open(pointdisp_filename, 'w') as f:
                    f.write(pointdisp_string)

    # noinspection PyMethodMayBeStatic
    def check_output_file(self, filename, nfaces):
        counter = 0
        nlines = 0
        lim = 1000000
        sleep_time = 0.01
        while (nlines < nfaces + 2) and counter < lim:
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

    def send_message(self, message):
        file = os.path.join(self.working_directory, message + '.coco')
        open(file, 'w').close()
        return

    def wait_message(self, message):
        wait_time_lim = 10 * 60  # 10 minutes maximum waiting time for a single flow solver iteration
        cumul_time = 0
        file = os.path.join(self.working_directory, message + '.coco')
        while not os.path.isfile(file):
            time.sleep(0.01)
            cumul_time += 0.01
            if cumul_time > wait_time_lim:
                self.openfoam_process.kill()
                self.openfoam_process.wait()
                raise RuntimeError(f'CoCoNuT timed out in the OpenFOAM solver_wrapper, waiting for message: '
                                   f'{message}.coco')
        os.remove(file)
        return

    def check_message(self, message):
        file = os.path.join(self.working_directory, message + '.coco')
        if os.path.isfile(file):
            os.remove(file)
            return True
        return False

    def remove_all_messages(self):
        for file_name in os.listdir(self.working_directory):
            if file_name.endswith('.coco'):
                file = os.path.join(self.working_directory, file_name)
                os.remove(file)

    def check_software(self):
        if subprocess.check_call(self.application + ' -help &> checkSoftware', shell=True, env=self.env) != 0:
            raise RuntimeError(
                f'OpenFOAM not loaded properly. Check if the solver load commands for the "machine_name" are correct.')

        # check version
        with open('checkSoftware', 'r') as f:
            last_line = f.readlines()[-2]  # second last line contains 'Build: XX' with XX the version number
        os.remove('checkSoftware')
        version_nr = last_line.split(' ')[-1]
        if version_nr[:-1] != self.version:
            raise RuntimeError(
                f'OpenFOAM-{self.version} should be loaded! Currently, OpenFOAM-{version_nr[:-1]} is loaded')

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
        self.write_interval = of_io.get_int(input_string=control_dict, keyword='writeInterval')
        time_format = of_io.get_string(input_string=control_dict, keyword='timeFormat')
        self.write_precision = of_io.get_int(input_string=control_dict, keyword='writePrecision')
        adjustTimeStep = of_io.get_string(input_string=control_dict, keyword='adjustTimeStep')
        if adjustTimeStep == "no":
            self.adjustFSITimeStep = False
        else:
            self.adjustFSITimeStep = True

        if not time_format == 'fixed':
            msg = f'timeFormat:{time_format} in controlDict not implemented. Changed to "fixed"'
            tools.print_info(msg, layout='warning')
            control_dict = re.sub(r'timeFormat' + of_io.delimter + r'\w+', f'timeFormat    fixed',
                                  control_dict)
        control_dict = re.sub(r'application' + of_io.delimter + r'\w+', f'{"application":<16}{self.application}',
                              control_dict)
        control_dict = re.sub(r'startTime' + of_io.delimter + of_io.float_pattern,
                              f'{"startTime":<16}{self.start_time}', control_dict)
        control_dict = re.sub(r'deltaT' + of_io.delimter + of_io.float_pattern, f'{"deltaT":<16}{self.delta_t}',
                              control_dict)
        control_dict = re.sub(r'timePrecision' + of_io.delimter + of_io.int_pattern,
                              f'{"timePrecision":<16}{self.time_precision}',
                              control_dict)
        control_dict = re.sub(r'endTime' + of_io.delimter + of_io.float_pattern, f'{"endTime":<16}1e15', control_dict)


        # delete previously defined coconut functions
        coconut_start_string = '\n// CoCoNuT function objects'
        control_dict = re.sub(coconut_start_string + r'.*', '', control_dict, flags=re.S)

        with open(file_name, 'w') as control_dict_file:
            control_dict_file.write(control_dict)
            control_dict_file.write(coconut_start_string + '\n')
            control_dict_file.write('boundary_names (')

            for boundary_name in self.boundary_names:
                control_dict_file.write(boundary_name + ' ')

            control_dict_file.write(');\n\n')
            control_dict_file.write('functions\n{\n')

            control_dict_file.write(self.wall_shear_stress_dict())
            for boundary_name in self.boundary_names:
                control_dict_file.write(self.pressure_and_traction_dict(boundary_name))
            control_dict_file.write('}')
            self.write_footer(file_name)

    def wall_shear_stress_dict(self):
        raise NotImplementedError('Base class method is called, should be implemented in derived class')

    def pressure_and_traction_dict(self, boundary_name):
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
                            raise RuntimeError(f'Variable: {variable} equation is not solved in {self.application}')

                    with open(self.res_filepath, 'a') as f:
                        np.savetxt(f, [residual_array], delimiter=', ')
