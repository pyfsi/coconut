from coconut import data_structure
from coconut.coupling_components.component import Component
from coconut.data_structure.interface import Interface
from coconut import tools
from subprocess import check_call

import numpy as np
import os
import linecache
import sys
import time
import copy
import subprocess


def create(parameters):
    return SolverWrapperOpenFOAM_41(parameters)


class SolverWrapperOpenFOAM_41(Component):
    def __init__(self, parameters):
        super().__init__()

        # Settings
        self.settings = parameters["settings"]
        self.working_directory = self.settings["working_directory"]
        self.moduleType = "OpenFOAM"  # hard-coded, cannot be altered by naive user
        self.moduleVersion = "4.1"  # hard-coded, cannot be altered by naive user
        self.module = self.moduleType + "/" + self.moduleVersion
        self.application = self.settings[
            "application"]  # What type of OF-solver to be used - solver requires adaptation before running with OpenFOAM - name of adapted solver starts with 'CoCoNuT_'
        self.dimensions = self.settings["dimensions"]
        if (self.dimensions != 2) and (self.dimensions != 3):
            sys.exit("OpenFOAM-case should be 2D or 3D.")
        self.dt = self.settings["dt"]  # Time step size
        self.start_time = self.settings[
            "start_time"]  # Start time - also the name of folder containing data at this time
        self.end_time = self.settings["end_time"]  # End time
        self.cores = self.settings["cores"]  # Number of cores to be used in the OpenFOAM-calculation
        self.decomposeMethod = self.settings[
            "decomposeMethod"]  # Decomposition-method, can be "simple", "scotch"
        self.newtonmax = self.settings["newtonmax"]  # Maximal number of Newton iterations
        self.newtontol = self.settings["newtontol"]  # Tolerance of Newton iterations
        self.write_interval = self.settings[
            "write_interval"]  # Number of time steps between consecutive saves performed by OpenFOAM
        self.write_precision = self.settings["write_precision"]  # writePrecision-parameter in OpenFOAM
        self.time_precision = self.settings["time_precision"]  # timePrecision-parameter in OpenFOAM
        self.boundary_names = self.settings[
            'boundary_names']  # boundary_names is the set of boundaries where the moving interface is located (will be used to go through OF-files)
        self.meshmotion_solver = self.settings["meshmotion_solver"]
        self.diffusivity = self.settings["diffusivity"]

        # debug
        self.debug = False  # set on True to save copy of input and output files in every iteration

        # Check that the boundary_names and the interface_input and interface_output are defined consistently in the JSON-file
        # For every boundary_name element, there should be one interface_input (boundary_name+"_input") element and one interface_output (boundary_name+"_output") element.
        # Make the distinction between both: interface_input/output are names of the pyKratos ModelPart - boundary_names is the name of the boundary as defined in OpenFOAM!
        self.check_interfaces()

        # Check that the correct modules have been loaded
        self.check_software()

        # Remove possible CoCoNuT-message from previous interrupt
        self.remove_all_messages()

        # Creating OpenFOAM-files - raw dictionary files are predefined in the solver_wrapper folder (and should not be moved)
        # DecomposeParDict: replace raw settings by actual settings defined by user in json-file
        if self.cores > 1:  # Only if calculating in parallel
            decomposeParDict_raw_name = os.path.join(os.path.realpath(os.path.dirname(__file__)),
                                                     "decomposeParDict_raw")
            decomposeParDict_name = os.path.join(self.working_directory, "system/decomposeParDict")
            with open(decomposeParDict_raw_name, 'r') as raw_file:
                with open(decomposeParDict_name, 'w') as new_file:
                    for line in raw_file:
                        line = line.replace('|CORES|', str(self.cores))
                        line = line.replace('|DECOMPOSEMETHOD|', str(self.decomposeMethod))
                        new_file.write(line)

            self.write_footer(decomposeParDict_name)

        # ControlDict: replace raw settings by actual settings defined by user in json-file AND add function objects to write pressure and wall shear stress
        controlDict_raw_name = os.path.join(os.path.realpath(os.path.dirname(__file__)), "controlDict_raw")
        controlDict_name = os.path.join(self.working_directory, "system/controlDict")
        with open(controlDict_raw_name, 'r') as raw_file:
            with open(controlDict_name, 'w') as new_file:
                for line in raw_file:
                    line = line.replace('|APPLICATION|', str(self.application))
                    line = line.replace('|START_TIME|', str(self.start_time))
                    line = line.replace('|END_TIME|', str(self.end_time))
                    line = line.replace('|DT|', str(self.dt))
                    line = line.replace('|WRITE_INTERVAL|', str(self.write_interval))
                    line = line.replace('|WRITE_PRECISION|', str(self.write_precision))
                    line = line.replace('|TIME_PRECISION|', str(self.time_precision))
                    if '|BOUNDARY_NAMES|' in line:
                        firstBoundary = True
                        for interfaces in self.boundary_names:
                            if firstBoundary:
                                boundary_name_temp = "(" + interfaces
                                firstBoundary = False
                            else:
                                boundary_name_temp += " , " + interfaces
                        boundary_name_temp += ")"
                        line = line.replace('|BOUNDARY_NAMES|', boundary_name_temp)
                    new_file.write(line)

        n_key = 0
        if len(self.boundary_names) == 1:
            for key in self.boundary_names:
                self.write_control_dict_function(controlDict_name, "surfaceRegion", "libfieldFunctionObjects",
                                                 "PRESSURE", key, True, False)
                self.write_control_dict_function(controlDict_name, "wallShearStress", "libfieldFunctionObjects",
                                                 "wallShearStress", key, False, False)
                self.write_control_dict_function(controlDict_name, "surfaceRegion", "libfieldFunctionObjects",
                                                 "TRACTION", key, False, True)
        else:
            for key in self.boundary_names:
                if n_key == 0:
                    self.write_control_dict_function(controlDict_name, "surfaceRegion", "libfieldFunctionObjects",
                                                     "PRESSURE", key, True, False)
                    self.write_control_dict_function(controlDict_name, "wallShearStress", "libfieldFunctionObjects",
                                                     "wallShearStress", key, False, False)
                    self.write_control_dict_function(controlDict_name, "surfaceRegion", "libfieldFunctionObjects",
                                                     "TRACTION", key, False, False)
                else:
                    self.write_control_dict_function(controlDict_name, "surfaceRegion", "libfieldFunctionObjects",
                                                     "PRESSURE", key, False, False)
                    if n_key == (len(self.boundary_names) - 1):
                        self.write_control_dict_function(controlDict_name, "surfaceRegion", "libfieldFunctionObjects",
                                                         "TRACTION", key, False, True)
                    else:
                        self.write_control_dict_function(controlDict_name, "surfaceRegion", "libfieldFunctionObjects",
                                                         "TRACTION", key, False, False)
                n_key += 1
        self.write_footer(controlDict_name)
        # DynamicMeshDict: replace raw settings by actual settings defined by user in json-file
        dynamicMeshDict_raw_name = os.path.join(os.path.realpath(os.path.dirname(__file__)), "dynamicMeshDict_raw")
        dynamic_mesh_dict_name = os.path.join(self.working_directory, "constant/dynamicMeshDict")
        str_boundary = ""
        for boundary in self.boundary_names:
            str_boundary = str_boundary + " " + str(boundary)
        with open(dynamicMeshDict_raw_name, 'r') as raw_file:
            with open(dynamic_mesh_dict_name, 'w') as new_file:
                for line in raw_file:
                    line = line.replace('|MESHMOTION_SOLVER|', str(self.meshmotion_solver))
                    line = line.replace('|DIFFUSIVITY|', str(self.diffusivity))
                    line = line.replace('|NUM_INTERFACE_INPUT|', str(len(self.settings['boundary_names'])))
                    line = line.replace('|INTERFACE_INPUT|', str_boundary)
                    new_file.write(line)

        self.write_footer(dynamic_mesh_dict_name)

        # Creating Model
        self.model = data_structure.Model()

        # Creating ModelParts and adding variables to these ModelParts - should happen before node addition
        # for key, value in (self.settings['interface_input'].items() + self.settings['interface_output'].items()):
        #     self.model.create_model_part(key)
        #     mp = self.model[key]
        #     for var_name in value:
        #         var = vars(data_structure)[var_name]
        #         mp.AddNodalSolutionStepVariable(var)

        # Adding nodes to ModelParts - should happen after variable definition; writeCellcentres writes cellcentres in internal field and face centres in boundaryField
        check_call("cd " + self.working_directory + "; writeCellCentres -time " + str(
            self.start_time) + " &> log.writeCellCentres;", shell=True)

        # Want "cellCentres for face boundaries"
        for boundary in self.boundary_names:
            source_file = self.working_directory + "/constant/polyMesh"
            node_ids, node_coords, face_centres, start_face, nfaces, self.total_nfaces = self.get_point_ids(boundary,
                                                                                                            source_file)
            # create input model part
            self.model.create_model_part(f'{boundary}_input', node_coords[:, 0], node_coords[:, 1], node_coords[:, 2],
                                         node_ids)

            name_x = os.path.join(self.working_directory, "0/ccx")
            name_y = os.path.join(self.working_directory, "0/ccy")
            name_z = os.path.join(self.working_directory, "0/ccz")
            index_x = self.find_string_in_file(boundary, name_x)
            index_y = self.find_string_in_file(boundary, name_y)
            index_z = self.find_string_in_file(boundary, name_z)

            with open(name_x, 'r') as fx:
                fx_lines = fx.readlines()

            with open(name_y, 'r') as fy:
                fy_lines = fy.readlines()

            with open(name_z, 'r') as fz:
                fz_lines = fz.readlines()

            x0 = np.empty(len(face_centres))
            y0 = np.empty(len(face_centres))
            z0 = np.empty(len(face_centres))
            ids = np.arange(0, len(face_centres))

            for i in np.arange(0, len(face_centres)):
                x0[i] = float(fx_lines[i + 6 + index_x].split("\n")[0])
                y0[i] = float(fy_lines[i + 6 + index_y].split("\n")[0])
                z0[i] = float(fz_lines[i + 6 + index_z].split("\n")[0])

            # create output model part
            mp_output = self.model.create_model_part(f'{boundary}_output', x0, y0, z0, ids)
            mp_output.start_face = start_face
            mp_output.nfaces = nfaces

        # Create CoSimulationInterfaces
        self.interface_input = Interface(self.settings["interface_input"], self.model)
        self.interface_output = Interface(self.settings["interface_output"], self.model)


    def initialize(self):
        super().initialize()

        # Define timestep and physical time
        self.timestep = 0
        self.physical_time = self.start_time

        # If no pointDisplacement file is defined yet, initialize a pointDisplacement file in the start time folder
        # Normally, after restart or from the second iteration onwards, a pointDisplacement-file already exists. In that case, that pointDisplacement-file will be used (and is NOT overwritten)
        pointdisp_raw_name = os.path.join(os.path.realpath(os.path.dirname(__file__)), "pointDisplacement_raw")
        pointdisp_name = os.path.join(self.working_directory, str(self.physical_time), "pointDisplacement")
        if not (os.path.isfile(pointdisp_name)):
            self.write_pointdisplacement_file(pointdisp_raw_name, pointdisp_name)

        ##Copy zero folder to folder with correctly named timeformat
        if self.physical_time == 0:
            timestamp = '{:.{}f}'.format(self.physical_time, self.time_precision)
            path_orig = os.path.join(self.working_directory, "0")
            path_new = os.path.join(self.working_directory, timestamp)
            os.system("cp -r " + path_orig + " " + path_new)

        ## If parallell do a decomposition and establish a remapping for the output based on the faceProcAddressing
        '''Note concerning the sequence: The file ./processorX/constant/polyMesh/pointprocAddressing contains a list of 
        indices referring to the original index in the ./constant/polyMesh/points file, these indices go from 0 to nPoints -1
        However, mesh faces can be shared between processors and it has to be tracekd whether these are inverted or not
        This inversion is indicated by negative indices. However, as minus 0 is not a thing, the indices are first incremented by 1 before inversion
        Therefore to get the correct index one should use |index|-1!!
        Normally no doubles should be encountered on an interface as these faces are not shared by processors
        '''
        if self.cores > 1:
            check_call("cd " + self.working_directory + "; decomposePar -force -time " + str(
                self.start_time) + " &> log.decomposePar;", shell=True)
            for boundary in self.boundary_names:
                mp_output = self.model.get_model_part(f'{boundary}_output')
                mp_output.sequence = []

                for p in range(self.cores):
                    count = 0
                    path = os.path.join(self.working_directory, "processor" + str(p), "constant", "polyMesh",
                                        "faceProcAddressing")
                    with open(path, 'r') as f:
                        face_Lines = f.readlines()
                    nfaces = int(face_Lines[18].split("\n")[0])
                    for i in range(20, 20 + nfaces):
                        ind = np.abs(int(face_Lines[i].split("\n")[0])) - 1

                        if ind >= mp_output.start_face and ind < mp_output.start_face + mp_output.nfaces:
                            mp_output.sequence.append(ind - mp_output.start_face)
                            count += 1

                np.savetxt(os.path.join(self.working_directory, f'sequence_{boundary}.txt'),
                           np.array(mp_output.sequence).astype(int), fmt='%i')
                if len(mp_output.sequence) != mp_output.nfaces:
                    print(f"sequence: {len(mp_output.sequence)}")
                    print(f"nNodes: {mp_output.NumberOfNodes}")
                    raise ValueError("Number of face indices in sequence does not correspond to number of elements")

        # Don't forget to start OpenFOAM-loop!
        if self.cores == 1:
            cmd = self.application + "&> log." + self.application
        else:
            cmd = "mpirun -np " + str(self.cores) + " " + self.application + " -parallel &> log." + self.application

        self.openfoam_process = subprocess.Popen(cmd, cwd=self.working_directory, shell=True)

        ### CoConuT_OpenFOAMSolver is running from here on

    def initialize_solution_step(self):
        super().initialize_solution_step()

        # For parallel runs need to create a folder with the correct time stamp for decomposition of pointDisplacement_Next
        # For serial runs, this folder will normally be present
        timestamp = '{:.{}f}'.format(self.physical_time, self.time_precision)
        path = os.path.join(self.working_directory, timestamp)
        if self.cores > 1:
            os.system('mkdir -p ' + path)
        elif self.physical_time == 0:  # for serial also need to make a folder 0.0000 with specified precision
            os.system('mkdir -p ' + path)

        # # The displacement of the FSI interface is passed through pointDisplacement_Next, which is prepared here
        # pointdisp_raw_name=os.path.join(os.path.realpath(os.path.dirname(__file__)),"pointDisplacement_raw")
        # pointdisp_name=os.path.join(self.working_directory,str(self.physical_time),timestamp,"pointDisplacement_Next")
        # self.write_pointdisplacement_file(pointdisp_raw_name,pointdisp_name)

        # Prepare new time step folder and reset the number of iterations
        self.timestep += 1
        self.iteration = 0
        self.physical_time += self.dt

        self.prev_timestamp = timestamp
        self.cur_timestamp = '{:.{}f}'.format(self.physical_time, self.time_precision)

        if self.cores < 2:  # if serial
            new_path = os.path.join(self.working_directory, self.cur_timestamp)
            if os.path.isdir(new_path):
                print("\n\n\n Warning! In 5s, CoCoNuT will overwrite existing time step folder: " + str(
                    new_path) + ". \n\n\n")
                time.sleep(5)
                os.system("rm -rf " + new_path)
            os.system("mkdir -p " + new_path)
        else:
            for i in np.arange(self.cores):
                new_path = os.path.join(self.working_directory, "processor" + str(i), self.cur_timestamp)
                if os.path.isdir(new_path):
                    if i == 0:
                        print(
                            "\n\n\n Warning! In 5s, CoCoNuT will overwrite existing time step folder in processor-subfolders. \n\n\n")
                        time.sleep(5)
                    os.system("rm -rf " + new_path)
                os.system("mkdir -p " + new_path)

        self.send_message('next')  # Let OpenFOAM go to next time step
        self.wait_message('next_ready')  # Let OpenFOAM wait for input data

    @tools.time_solve_solution_step
    def solve_solution_step(self, interface_input):  # NOT CHANGED YET! PURELY COPIED FROM FLUENT WRAPPER!!!!!!
        self.iteration += 1

        # store incoming displacements
        self.interface_input.set_interface_data(interface_input.get_interface_data())

        # write interface data to OpenFOAM-file
        self.write_node_input()

        # copy output data for debugging
        if self.debug:
            if self.cores > 1:
                for i in range(0, self.cores):
                    path = os.path.join(self.working_directory, "processor" + str(i), self.prev_timestamp,
                                        "pointDisplacement_Next")
                    path2 = os.path.join(self.working_directory, "processor" + str(i), self.prev_timestamp,
                                         "pointDisplacement_Next_Iter" + str(self.iteration))
                    cmd = f"cp {path} {path2}"
                    os.system(cmd)
            else:
                path = os.path.join(self.working_directory, self.prev_timestamp, "pointDisplacement_Next")
                path2 = os.path.join(self.working_directory, self.prev_timestamp,
                                     "pointDisplacement_Next_Iter" + str(self.iteration))
                cmd = f"cp {path} {path2}"
                os.system(cmd)

        # let OpenFOAM run, wait for data
        '''OpenFOAM tends to keep on appending to files while already providing access, this causes issues when reading out the data
        Therefore these files are removed if they already exist prior to generating new data output '''
        for boundary in self.boundary_names:
            # specify location of pressure and traction
            traction_name = "TRACTION_" + boundary
            pressure_name = "PRESSURE_" + boundary
            wss_file = os.path.join(self.working_directory, "postProcessing", traction_name, "surface",
                                    self.cur_timestamp, "wallShearStress_patch_" + boundary + ".raw")
            pres_file = os.path.join(self.working_directory, "postProcessing", pressure_name, "surface",
                                     self.cur_timestamp, "p_patch_" + boundary + ".raw")
            if os.path.isfile(wss_file):
                os.system(f"rm -rf {wss_file}")
            if os.path.isfile(pres_file):
                os.system(f"rm -rf {pres_file}")

        self.send_message('continue')
        self.wait_message('continue_ready')

        # read data from OpenFOAM
        self.read_node_output()

        # return interface_output object
        return self.interface_output

    def finalize_solution_step(self):
        super().finalize_solution_step()

        # Let OpenFOAM check whether it needs to save this timestep (in OF-solver: runTime.write())

        if not (self.timestep % self.write_interval):
            if self.cores > 1:  # Remove folder that was used for pointDisplacement_Next
                # at end of time step if parallel run if not writeInterval
                path = os.path.join(self.working_directory, self.prev_timestamp)
                os.system("rm -rf " + path)
            self.send_message('save')
            self.wait_message('save_ready')

    def finalize(self):
        super().finalize()

        self.send_message('stop')
        self.wait_message('stop_ready')
        os.system(f"pkill -f {self.application}")
        self.openfoam_process.kill()


    def get_interface_input(self):
        return self.interface_input

    def get_interface_output(self):
        return self.interface_output

    def write_header(self, file_loc, class_name, object_name):
        f = open(file_loc, 'w')
        f.write(r'/*--------------------------------*- C++ -*----------------------------------*\\' + "\n")
        f.write(r'| =========                 |                                                 |' + "\n")
        f.write(r'| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |' + "\n")
        f.write(r'|  \\    /   O peration     | Version:  4.x                                   |' + "\n")
        f.write(r'|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |' + "\n")
        f.write(r'|    \\/     M anipulation  |                                                 |' + "\n")
        f.write(r'\*---------------------------------------------------------------------------*/' + "\n")
        f.write(r'FoamFile' + "\n")
        f.write(r'{' + "\n")
        f.write('\t version \t\t 4.1;' + "\n")
        f.write('\t format \t\t ascii;' + "\n")
        f.write('\t class \t\t ' + class_name + ';' + "\n")
        f.write('\t object \t\t ' + object_name + ';' + "\n")
        f.write('}' + "\n")
        f.write(r'// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //' + "\n")
        f.write("\n")
        f.close()

    def read_node_output(self):
        ''' This is to be verified but it might be imortant that when a calculation is started from a prior time step while there is still
        data in the postprocessing folder for that time step from a previous run, that then the checks to see whether the file is completely updated might fail
        and consequently cause the simulation to crash
        Maybe it would be better to remove these files beforehand if they exist instead of overwriting them
        This however can not be done here, but should be done in another part of the code '''
        n_key = 0

        # Default value is 1.0 for compressible case. When the solver is incompressible, the pressure and shear stress are
        # kinematic; therefore multiply with the fluid density.
        density = 1.0
        if self.settings['is_incompressible']:
            density = self.settings['density']

        for boundary in self.boundary_names:
            # specify location of pressure and traction
            traction_name = "TRACTION_" + boundary
            pressure_name = "PRESSURE_" + boundary
            mp_name = f'{boundary}_output'
            mp = self.model.get_model_part(mp_name)
            nfaces_tot = mp.size
            wss_tmp = np.zeros([nfaces_tot, 3])
            pres_tmp = np.zeros([nfaces_tot, 1])
            wss_file = os.path.join(self.working_directory, "postProcessing", traction_name, "surface",
                                    self.cur_timestamp, "wallShearStress_patch_" + boundary + ".raw")
            pres_file = os.path.join(self.working_directory, "postProcessing", pressure_name, "surface",
                                     self.cur_timestamp, "p_patch_" + boundary + ".raw")

            # read traction
            counter = 0
            nlines = 0
            lim = 1000
            while (nlines < nfaces_tot + 2) and counter < lim:
                if os.path.isfile(wss_file):
                    nlines = sum(1 for line in open(wss_file))
                time.sleep(0.01)
                counter += 1
            if counter == lim:
                raise RuntimeError("Timed out waiting for wss file: " + wss_file)

            with open(wss_file, 'r') as f:
                f_lines = f.readlines()
            index_start = 2
            for i in range(nfaces_tot):
                split = f_lines[index_start + i].split()
                if self.physical_time == self.start_time:  # At start perform check of the mapping, CARE for restarts this check will need to be updated as the file will give coordinates of the deformed geometry
                    # Could use the displacement vector for this which we most likely need to read anyways
                    pos = np.array([float(split[0]), float(split[1]), float(split[2])])
                    pos_mp = np.array([mp.x0[mp.sequence[i]], mp.y0[mp.sequence[i]], mp.z0[mp.sequence[i]]])

                    if (np.linalg.norm(pos_mp - pos) > 1e-05):
                        raise ValueError("Positions do not agree !!")
                wss_tmp[i, 0] = split[3]
                wss_tmp[i, 1] = split[4]
                wss_tmp[i, 2] = split[5].split("\n")[0]
            f.close()

            # read pressure
            counter = 0
            nlines = 0
            while (nlines < nfaces_tot + 2) and counter < lim:
                nlines = 0
                if os.path.isfile(pres_file):
                    nlines = sum(1 for line in open(pres_file))

                time.sleep(0.01)
                counter += 1

            if counter == lim:
                raise RuntimeError("Timed out waiting for pressure file: " + pres_file)

            with open(pres_file, 'r') as f:
                f_lines = f.readlines()
            index_start = 2

            for i in np.arange(nfaces_tot):
                val = f_lines[index_start + i].split()[3].split("\n")[0]
                pres_tmp[i, 0] = float(val)
            f.close()

            if self.cores > 1:
                pos_list = mp.sequence
            else:
                pos_list = [pos for pos in range(0, nfaces_tot)]

            wall_shear_stress = np.empty(wss_tmp.shape)
            pressure = np.empty(pres_tmp.shape)

            wall_shear_stress[pos_list,] = wss_tmp[:, ]
            pressure[pos_list] = pres_tmp[:]

            self.interface_output.set_variable_data(mp_name, 'traction', wall_shear_stress * -1 * density)
            self.interface_output.set_variable_data(mp_name, 'pressure', pressure * density)

            # go to next interface
            n_key += 1

    # writeFooter: to write OpenFOAM-footer in file at location 'file_loc'
    def write_footer(self, file_loc):
        f = open(file_loc, 'a+')
        f.write("\n")
        f.write(r'// ************************************************************************* //' + "\n")
        f.close()

    def write_node_input(self):
        # The displacement of the FSI interface is passed through the file "pointDisplacement_Next"
        # This function will prepare that file in a "serial format" and then decompose it for parallel operation

        pointdisp_raw_name = os.path.join(os.path.realpath(os.path.dirname(__file__)), "pointDisplacement_raw")
        pointdisp_name = os.path.join(self.working_directory, self.prev_timestamp, "pointDisplacement_Next")
        self.write_pointdisplacement_file(pointdisp_raw_name, pointdisp_name)

        disp_file = pointdisp_name

        n_key = 0
        for boundary in self.boundary_names:
            mp_name = f'{boundary}_input'
            displacement = self.interface_input.get_variable_data(mp_name, 'displacement')

            start_nr = self.find_string_in_file(boundary, disp_file)
            os.system("head -n " + str(start_nr + 1) + " " + disp_file + " > tempDisp")

            with open('tempDisp', 'a+') as file:
                file.write("\t { \n")
                file.write("\t\t type  \t fixedValue; \n")
                file.write('\t\t value \t nonuniform List<vector> ( \n')
                for i in range(displacement.shape[0]):
                    file.write(' (' + f'{displacement[i, 0]:27.17e} {displacement[i, 1]:27.17e} {displacement[
                        i, 2]:27.17e}' + ') \n')
                file.write(');\n')

            os.system("wc -l " + disp_file + " > lengthDisp")
            lengthDisp_file = open("lengthDisp", 'r')
            length_disp = int(lengthDisp_file.readline().split(" ")[0])
            lengthDisp_file.close()
            os.system("tail -n " + str(length_disp - (start_nr + 1)) + " " + disp_file + " > tempDisp2")
            start_to_end_nr = self.find_string_in_file("}", "tempDisp2")
            os.system(
                "tail -n " + str(length_disp - (start_nr + 1) - start_to_end_nr) + " " + disp_file + " > tempDisp3")
            os.system("cat tempDisp tempDisp3 > " + disp_file)
            os.system("rm tempDisp* lengthDisp")
            n_key += 1

        if self.cores > 1:
            check_call(
                "cd " + self.working_directory + "; decomposePar -fields -time " + self.prev_timestamp + " &> log.decomposePar;",
                shell=True)

    def write_control_dict_function(self, file_name, func_name, lib_name, variable_name, patch_name, write_start,
                                    write_end):
        with open(file_name, 'a+') as file:
            if write_start:
                file.write("functions \n")
                file.write("{ \n ")
            if variable_name == "wallShearStress":
                file.write(" \n \t " + variable_name + " \n")
            else:
                file.write(" \n \t " + variable_name + "_" + patch_name + " \n")
            file.write("\t { \n")
            file.write("\t\t type  \t " + func_name + "; \n")
            file.write('\t\t libs \t ("' + lib_name + '.so"); \n')
            file.write('\t\t executeControl \t timeStep; \n')
            file.write('\t\t executeInterval \t 1; \n')
            file.write('\t\t writeControl \t timeStep; \n')
            file.write('\t\t writeInterval \t 1; \n')
            file.write('\t\t timeFormat \t fixed; \n')
            file.write(f'\t\t timePrecision \t {self.time_precision}; \n')
            if func_name == "surfaceRegion":
                file.write('\t\t operation \t none; \n')
                file.write('\t\t writeFields \t true; \n')
                file.write('\t\t surfaceFormat \t raw; \n')
                file.write('\t\t regionType \t patch; \n')
                file.write('\t\t name \t ' + patch_name + ' ; \n')
                file.write('\t\t fields \n')
                file.write('\t\t ( \n')
                if variable_name == "PRESSURE":
                    file.write('\t\t\t p \n ')
                elif variable_name == "TRACTION":
                    file.write('\t\t\t wallShearStress \n')
                file.write("\t\t ); \n")
            elif func_name == "wallShearStress":
                file.write('\t\t log \t false; \n')
            file.write("\t } \n\n")
            if write_end:
                file.write("} \n ")

        file.close()

    def write_pointdisplacement_file(self, pointdisp_raw_name, pointdisp_name):
        with open(pointdisp_raw_name, 'r') as raw_file:
            with open(pointdisp_name, 'w') as new_file:
                for line in raw_file:
                    new_file.write(line)

                for boundary in self.boundary_names:
                    new_file.write(" \n \t " + boundary + " \n")
                    new_file.write("\t { \n")
                    new_file.write("\t\t type  \t fixedValue; \n")
                    new_file.write("\t\t value \t uniform (0 0 0); \n")
                    new_file.write("\t } \n")
                new_file.write("} \n")
                new_file.close()
                self.write_footer(pointdisp_name)
        raw_file.close()

    def find_string_in_file(self, string_name, file_name):
        index = -1
        with open(file_name) as f:
            for num, line in enumerate(f):
                if string_name in line:
                    index = num
                    break
        f.close()
        return index

    def send_message(self, message):
        file = os.path.join(self.working_directory, message + ".coco")
        open(file, 'w').close()
        return

    def wait_message(self, message):
        waitTimeLimit = 10 * 60  # 10 minutes maximum waiting time for a single flow solver iteration
        cumulTime = 0
        file = os.path.join(self.working_directory, message + ".coco")
        while not os.path.isfile(file):
            time.sleep(0.01)
            cumulTime += 0.01
            if cumulTime > waitTimeLimit:
                os.system("pkill " + self.application)
                sys.exit("CoCoNuT timed out in the OpenFOAM solver_wrapper, waiting for message: " + message + ".coco.")
        os.remove(file)
        return

    def check_message(self, message):
        file = os.path.join(self.working_directory, message + ".coco")
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
        if os.system(self.application + ' -help &> checkSoftware') != 0:
            sys.exit(
                "You either did not load the module for OpenFOAM/4.1, did not compile the solver you are trying to use or did not source $FOAM_BASH prior to execution of CoCoNuT.")

        # The statement above could work for other version of OpenFOAM is the solver is also compiled for that version. Therefore, the version number is checked explicity (don't forget to remove the end-of-line variable at the end of the String versionNr
        with open('checkSoftware', 'r') as f:
            lastLine = f.readlines()[-2]  # Second last line contains 'Build: XX' with XX the version number
        f.close()
        os.system('rm checkSoftware')
        versionNr = lastLine.split(' ')[-1]
        if versionNr[:-1] != self.moduleVersion:
            sys.exit("OpenFOAM 4.1 should be loaded! Currently, another version of OpenFOAM is loaded")

    def get_point_ids(self, boundary, dir):
        'Function that returns the local point IDs belonging to a specified boundary in the correct sequence for the pointDisplacement file'
        f_b = f'{dir}/boundary'
        f_f = f'{dir}/faces'
        f_p = f'{dir}/points'

        # Identify startface and endface for the boundary
        with open(f_b) as f:
            line = f.readline()
            while not boundary in line:
                line = f.readline()
            for i in range(0, 4):
                line = f.readline()
            nfaces = int(line[:-2].split()[1])
            line = f.readline()
            start_face = int(line[:-2].split()[1])

        # Get number of points to keep a list of booleans
        prev_line = "("
        with open(f_p) as f:
            line = f.readline()
            while not "(" in line:
                prev_line = line
                line = f.readline()
            n_points = int(prev_line)
            points = np.zeros((n_points, 3))
            count = 0
            line = f.readline()
            while any(char.isdigit() for char in line):
                temp = line.split(" ")
                points[count, 0] = float(temp[0][1:])
                points[count, 1] = float(temp[1][:])
                points[count, 2] = float(temp[2][:-2])
                count += 1
                line = f.readline()

        points_bool = np.zeros((n_points, 1))
        boundary_ind = []

        # Read in nodes file

        # Read in the list of faces and the nodes constituting those faces
        all_face_nodes = []
        with open(f_f) as f:
            line = f.readline()
            while not "(" in line:
                line = f.readline()
            line = f.readline()
            while any(char.isdigit() for char in line):
                list = line[2:-2].split()
                all_face_nodes.append(list)
                line = f.readline()

        # Extract the cell faces belonging to the boundary and create an ordered list of node id's
        face_centers = np.zeros((nfaces, 3))
        node_coords = []
        for nf in range(start_face, start_face + nfaces):
            face_center = np.zeros([1, 3])
            list = all_face_nodes[nf]
            for n in list:
                face_center += points[int(n), :]
                if points_bool[int(n) - 1, 0] == 0:  # If not yet in list, add it to the list
                    boundary_ind.append(int(n))
                    node_coords.append(points[int(n), :])
                    points_bool[int(n) - 1, 0] = 1
            face_centers[nf - start_face, :] = face_center / float(len(list))
        boundary_ind = np.array(boundary_ind)
        node_coords = np.array(node_coords, dtype=float)
        os.path.join(self.working_directory, "faceCenters.txt")
        np.savetxt(os.path.join(self.working_directory, "faceCenters.txt"), face_centers)

        return boundary_ind, node_coords, face_centers, start_face, nfaces, len(all_face_nodes)

    def check_interfaces(self):
        input_interface_model_parts = [param["model_part"] for param in self.settings["interface_input"]]
        output_interface_model_parts = [param["model_part"] for param in self.settings["interface_output"]]
        boundary_names = self.settings["boundary_names"]

        for boundary_name in boundary_names:
            if not f'{boundary_name}_input' in input_interface_model_parts:
                raise RuntimeError(
                    f'Error in json file: {boundary_name}_input not listed in "interface_input": {self.settings[
                        "interface_input"]}.\n. <boundary> in the "boundary_names" in json file should have corresponding <boundary>_input in "interface_input" list. ')

            if not f'{boundary_name}_output' in output_interface_model_parts:
                raise RuntimeError(
                    f'Error in json file: {boundary_name}_output not listed in "interface_output": {self.settings[
                        "interface_output"]}.\n. <boundary> in the "boundary_names" in json file should have corresponding <boundary>_input in "interface_output" list.')
