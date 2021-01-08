from coconut import data_structure
from coconut.coupling_components.component import Component
from coconut.coupling_components import tools

import os
from os.path import join
import subprocess
import time
import numpy as np
import re


def create(parameters):
    return SolverWrapperAbaqus614(parameters)


class SolverWrapperAbaqus614(Component):
    def __init__(self, parameters):
        super().__init__()
        # settings
        """
               settings of solver_wrappers.abaqus.614:

                    working_directory       Absolute path to working directory
                                            or relative path w.r.t current directory
                    cores                   Number of cpus to be used by Abaqus
                    input_file              Name of the Abaqus input file (located in the directory where Python is 
                                            launched)
                    dimensions              dimensionality of the problem (2 or 3)
                    arraysize               declare a sufficiently large array size for load array in FORTRAN
                    CSM_dir                 relative path to directory for the files and execution of the structural 
                                            solver 
                    ramp                    Boolean: 0 for step load, 1 for ramp load in Abaqus
                    delta_t                 Time step size
                    timestep_start          Time step from which is to be started (initial = 0)
                    surface_IDS             List of the names of the surfaces that take part in the FSI, as they are 
                                            known by Abaqus
                    interface_input         Interface for the load points and their corresponding variables (pressure,
                                            traction).
                    interface_output        Interface for the output nodes and their corresponding variable(s)
                                            displacements)
                    mp_mode                 Mode of the parallel computing (currently only THREADS is accepted)
                    save_iterations         number of time steps between consecutive saves of Abaqus-files
        """

        self.settings = parameters['settings']
        self.dir_csm = join(os.getcwd(), self.settings['working_directory'])  # *** alternative for getcwd?
        path_src = os.path.realpath(os.path.dirname(__file__))

        self.cores = self.settings['cores']  # number of CPUs Abaqus has to use
        self.dimensions = self.settings['dimensions']
        if self.dimensions == 2:
            tools.print_info("Warning for axisymmetric cases:\n\tIn Abaqus these have to be constructed around the "
                             "y-axis. \n\tSwitching of x and y-coordinates might be necessary but should be "
                             "accomplished by using an appropriate mapper.", layout='warning')
        self.array_size = self.settings["arraysize"]
        self.delta_t = self.settings["delta_t"]
        self.timestep_start = self.settings["timestep_start"]
        self.surfaceIDs = self.settings['surfaceIDs']
        self.n_surfaces = len(self.surfaceIDs)
        self.thread_ids = [i for i in range(0, self.n_surfaces)]  # list(range(self.n_surfaces))?
        self.mp_mode = self.settings["mp_mode"]
        self.input_file = self.settings["input_file"]
        self.timestep = self.timestep_start
        self.iteration = None
        self.model_part_surface_ids = {}  # surface IDs corresponding to ModelParts

        if "subcycling" in self.settings.keys():
            self.subcycling = self.settings["subcycling"]
            if self.subcycling:
                self.minInc = self.settings["minInc"]
                self.initialInc = self.settings["initialInc"]
                self.maxNumInc = self.settings["maxNumInc"]
                self.maxInc = self.settings["maxInc"]
                self.ramp = self.settings["ramp"]
            else:
                self.maxNumInc = 1
                self.maxInc = self.delta_t
                self.ramp = 0
        else:
            self.subcycling = 0
            self.maxNumInc = 1
            self.ramp = 0

        # Upon (re)starting Abaqus needs to run USRInit.f
        # A restart requires Abaqus to be booted with a restart file

        # prepare abaqus_v6.env
        self.hostnames = []
        self.hostnames_unique = []
        with open(join(self.dir_csm, "AbaqusHosts.txt"), "r") as hostfile:
            for line in hostfile:
                self.hostnames.append(line.rstrip())
                if not line.rstrip() in self.hostnames_unique:
                    self.hostnames_unique.append(line.rstrip())
        self.hostname_replace = ""
        for hostname in self.hostnames_unique:
            self.hostname_replace += "[\'" + hostname + "\', " + str(self.hostnames.count(hostname)) + "], "
        self.hostname_replace = self.hostname_replace.rstrip(", ")
        with open(join(path_src, "abaqus_v6.env"), "r") as infile:
            with open(join(self.dir_csm, "abaqus_v6.env"), "w") as outfile:
                for line in infile:
                    line = line.replace("|HOSTNAME|", self.hostname_replace)
                    line = line.replace("|MP_MODE|", self.mp_mode)
                    line = line.replace("|PID|", str(os.getpid()))
                    line = line.replace("|PWD|", os.path.abspath(os.path.join(self.dir_csm, os.pardir)))
                    line = line.replace("|CSM_dir|", self.settings["working_directory"])
                    if "|" in line:
                        raise ValueError(f"The following line in abaqus_v6.env still contains a \"|\" after "
                                         f"substitution: \n \t{line} \n Probably a parameter was not substituted")
                    outfile.write(line)

        # Create start and restart file
        self.write_start_and_restart_inp(join(self.dir_csm, self.input_file), self.dir_csm + "/CSM_Time0.inp",
                                         self.dir_csm + "/CSM_Restart.inp")

        # Prepare Abaqus USRInit.f
        usr = "USRInit.f"
        with open(join(path_src, usr), "r") as infile:
            with open(join(self.dir_csm, "usrInit.f"), "w") as outfile:
                for line in infile:
                    line = line.replace("|dimension|", str(self.dimensions))
                    line = line.replace("|surfaces|", str(self.n_surfaces))
                    line = line.replace("|cpus|", str(self.cores))

                    # if PWD is too long then FORTRAN code can not compile so this needs special treatment
                    line = self.replace_fortran(line, "|PWD|", os.path.abspath(os.getcwd()))
                    line = self.replace_fortran(line, "|CSM_dir|", self.settings["working_directory"])
                    if "|" in line:
                        raise ValueError(f"The following line in USRInit.f still contains a \"|\" after substitution: "
                                         f"\n \t{line} \n Probably a parameter was not substituted")
                    outfile.write(line)

        # Compile Abaqus USRInit.f
        path_libusr = join(self.dir_csm, "libusr/")
        os.system("rm -rf " + path_libusr)
        os.system("mkdir " + path_libusr)
        cmd = "abaqus make library=usrInit.f directory=" + path_libusr + " >> AbaqusSolver.log 2>&1"
        commands = [cmd]
        self.run_shell(self.dir_csm, commands, name='Compile_USRInit')

        # Get load points from usrInit.f
        if self.timestep_start == 0:
            cmd1 = f"export PBS_NODEFILE=AbaqusHosts.txt && unset SLURM_GTIDS"  # To get this to work on HPC?
            cmd2 = f"rm -f CSM_Time{self.timestep_start}Surface*Faces.dat " \
                f"CSM_Time{self.timestep_start}Surface*FacesBis.dat"
            # The output files will have a name with a higher time step  ("job=") than the input file ("input=")
            cmd3 = f"abaqus job=CSM_Time{self.timestep_start + 1} input=CSM_Time{self.timestep_start} " \
                f"cpus=1 user=usrInit.f output_precision=full interactive >> AbaqusSolver.log 2>&1"
            commands = [cmd1, cmd2, cmd3]
            self.run_shell(self.dir_csm, commands, name='Abaqus_USRInit_Time0')
        else:
            cmd1 = f"export PBS_NODEFILE=AbaqusHosts.txt && unset SLURM_GTIDS"  # To get this to work on HPC?
            cmd2 = f"rm -f CSM_Time{self.timestep_start}Surface*Faces.dat " \
                f"CSM_Time{self.timestep_start}Surface*FacesBis.dat"
            cmd3 = f"abaqus job=CSM_Time{self.timestep_start + 1} oldjob=CSM_Time{self.timestep_start} " \
                f"input=CSM_Restart cpus=1 user=usrInit.f output_precision=full interactive >> AbaqusSolver.log 2>&1"
            commands = [cmd1, cmd2, cmd3]
            self.run_shell(self.dir_csm, commands, name=f'Abaqus_USRInit_Restart')

        # prepare GetOutput.cpp
        get_output = "GetOutput.cpp"
        temp_str = ""
        for j in range(0, self.n_surfaces - 1):
            temp_str += f"\"{self.surfaceIDs[j]}\", "
        temp_str += f"\"{self.surfaceIDs[self.n_surfaces-1]}\""

        with open(join(path_src, get_output), "r") as infile:
            with open(join(self.dir_csm, get_output), "w") as outfile:
                for line in infile:
                    line = line.replace("|surfaces|", str(self.n_surfaces))
                    line = line.replace("|surfaceIDs|", temp_str)
                    line = line.replace("|dimension|", str(self.dimensions))
                    if "|" in line:
                        raise ValueError(f"The following line in GetOutput.cpp still contains a \"|\" after "
                                         f"substitution: \n \t{line} \n Probably a parameter was not subsituted")
                    outfile.write(line)

        # compile GetOutput.cpp
        cmd = "abaqus make job=GetOutput user=GetOutput.cpp >> AbaqusSolver.log 2>&1"
        commands = [cmd]
        self.run_shell(self.dir_csm, commands, name='Compile_GetOutput')

        # Get node positions (not load points) at timestep_start (0 is an argument to GetOutput.exe)
        cmd = f"abaqus ./GetOutput.exe CSM_Time{self.timestep_start+1} 0 >> AbaqusSolver.log 2>&1"
        commands = [cmd]
        self.run_shell(self.dir_csm, commands, name='GetOutput_Start')

        for i in range(0, self.n_surfaces):
            path_output = f"CSM_Time{self.timestep_start+1}Surface{i}Output.dat"
            path_nodes = f"CSM_Time{self.timestep_start}Surface{i}Nodes.dat"
            cmd = f"mv {path_output} {path_nodes}"
            commands = [cmd]
            self.run_shell(self.dir_csm, commands, name='Move_Output_File_To_Node')

            # Create elements file per surface
            face_file = os.path.join(self.dir_csm, f"CSM_Time{self.timestep_start}Surface{i}Cpu0Faces.dat")
            output_file = os.path.join(self.dir_csm, f"CSM_Time{self.timestep_start}Surface{i}Elements.dat")
            self.make_elements(face_file, output_file)

        # prepare Abaqus USR.f
        usr = "USR.f"
        with open(join(path_src, usr), "r") as infile:
            with open(join(self.dir_csm, "usr.f"), "w") as outfile:
                for line in infile:
                    line = line.replace("|dimension|", str(self.dimensions))
                    line = line.replace("|arraySize|", str(self.array_size))
                    line = line.replace("|surfaces|", str(self.n_surfaces))
                    line = line.replace("|cpus|", str(self.cores))
                    line = line.replace("|ramp|", str(self.ramp))
                    line = line.replace("|deltaT|", str(self.delta_t))

                    # if PWD is too long then FORTRAN code can not compile so this needs special treatment
                    line = self.replace_fortran(line, "|PWD|", os.path.abspath(os.getcwd()))
                    line = self.replace_fortran(line, "|CSM_dir|", self.settings["working_directory"])
                    if "|" in line:
                        raise ValueError(f"The following line in USR.f still contains a \"|\" after substitution: "
                                         f"\n \t{line} \n Probably a parameter was not subsituted")
                    outfile.write(line)

        # compile Abaqus USR.f
        os.system("rm -r " + path_libusr)  # remove libusr containing compiled USRInit.f
        os.system("mkdir " + path_libusr)
        cmd = "abaqus make library=usr.f directory=" + path_libusr + " >> AbaqusSolver.log 2>&1"
        commands = [cmd]
        self.run_shell(self.dir_csm, commands, name='Compile_USR')

        # create Model
        self.model = data_structure.Model()

        # create input ModelParts (load points)
        for item in (self.settings['interface_input']):
            mp_name = item['model_part']

            for i, surfaceID in enumerate(self.surfaceIDs):
                if surfaceID in mp_name:
                    self.model_part_surface_ids[mp_name] = i
                    break
            if mp_name not in self.model_part_surface_ids:
                raise AttributeError(f'Could not identify surfaceID corresponding to ModelPart {mp_name}.')
            mp_id = self.model_part_surface_ids[mp_name]

            # read in elements file
            tmp = f'CSM_Time{self.timestep_start}Surface{mp_id}Elements.dat'
            elem_file = join(self.dir_csm, tmp)
            elements = np.loadtxt(elem_file)

            if self.timestep_start > 0:
                tmp0 = f'CSM_Time0Surface{mp_id}Elements.dat'
                elem0_file = join(self.dir_csm, tmp0)
                elements0 = np.loadtxt(elem0_file)
            elif self.timestep_start == 0:
                elements0 = elements

            n_elem = int(elements[0])  # elements line 1 contains number of elements
            n_lp = int(elements[1])  # elements line 2 contains number of load points per element
            if elements.shape[0]-2 != int(n_elem):  # elements remainder contains element numbers involved in interface
                raise ValueError(f"Number of lines ({elements.shape[0]}) in {elem_file} does not correspond with the "
                                 f"number of elements ({n_elem})")

            if int(elements0[0]) != n_elem or int(elements0[1]) != n_lp:
                raise ValueError(f"Number of load points has changed for {mp_name}")

            # read in Faces file for load points and also file at time 0 for original positions for the mappers
            # TODO: do we need the current list of faces?
            tmp = f'CSM_Time{self.timestep_start}Surface{mp_id}Cpu0Faces.dat'
            faces_file = join(self.dir_csm, tmp)
            faces = np.loadtxt(faces_file)

            if self.timestep_start > 0:
                tmp0 = f'CSM_Time0Surface{mp_id}Cpu0Faces.dat'
                faces0_file = join(self.dir_csm, tmp0)
                faces0 = np.loadtxt(faces0_file)
            elif self.timestep_start == 0:
                faces0 = faces

            if faces.shape[1] != self.dimensions + 2:
                raise ValueError(f'given dimension does not match coordinates')

            # get load point coordinates and ids of load points
            prev_elem = 0
            prev_lp = 0
            ids = np.arange(n_elem*n_lp)
            coords_tmp = np.zeros((n_elem * n_lp, 3)).astype(float)  # z-coordinate mandatory: 0.0 for 2D
            for i in range(0, n_elem*n_lp):
                elem = int(faces0[i, 0])
                lp = int(faces0[i, 1])
                if elem < prev_elem:
                    raise ValueError(f"Element sequence is wrong ({elem}<{prev_elem})")
                elif elem == prev_elem and lp != prev_lp+1:
                    raise ValueError(f"Next line for same element ({elem}) does not contain next load point")
                elif elem > prev_elem and lp != 1:
                    raise ValueError(f"First line for element ({elem}) does not contain its first load point")
                if lp > n_lp:
                    raise ValueError(f"Load point ({lp}) exceeds the number of load points per element {n_lp}")

                # ids_tmp[i] = f"{elem}_{lp}"
                coords_tmp[i, :self.dimensions] = faces0[i, -self.dimensions:]  # extract last "dimensions" columns

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

            for i, surfaceID in enumerate(self.surfaceIDs):
                if surfaceID in mp_name:
                    self.model_part_surface_ids[mp_name] = i
                    break
            if mp_name not in self.model_part_surface_ids:
                raise AttributeError(f'Could not identify surfaceID corresponding to ModelPart {mp_name}.')
            mp_id = self.model_part_surface_ids[mp_name]

            # read in Nodes0 file
            tmp = f'CSM_Time{self.timestep_start}Surface{mp_id}Nodes.dat'
            nodes_file = join(self.dir_csm, tmp)
            nodes = np.loadtxt(nodes_file, skiprows=1)  # first line is a header

            if self.timestep_start > 0:
                tmp0 = f'CSM_Time0Surface{mp_id}Nodes.dat'
                nodes0_file = join(self.dir_csm, tmp0)
                nodes0 = np.loadtxt(nodes0_file, skiprows=1)
            elif self.timestep_start == 0:
                nodes0 = nodes

            if nodes.shape[1] != self.dimensions:
                raise ValueError(f'Given dimension does not match coordinates.')

            # get geometrical node coordinates
            n_nodes = nodes.shape[0]
            n_nodes0 = nodes0.shape[0]
            if n_nodes != n_nodes0:
                raise ValueError(f"Number of interface nodes has changed for {mp_name}")

            ids = np.arange(n_nodes)  # Abaqus does not use node ids but maintains the output order
            coords_tmp = np.zeros((n_nodes, 3)).astype(float)  # z-coordinate mandatory: 0.0 for 2D
            coords_tmp[:, :self.dimensions] = nodes0

            x0 = coords_tmp[:, 0]
            y0 = coords_tmp[:, 1]
            z0 = coords_tmp[:, 2]

            # create ModelPart
            self.model.create_model_part(mp_name, x0, y0, z0, ids)

        # check whether the ModelParts and output ModelParts have proper overlap
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
                tmp = f"CSM_Time{self.timestep}Surface{mp_id}Cpu0Input.dat"
                tmp2 = f"CSM_Time{self.timestep}Surface{mp_id}Cpu0Input_Iter{self.iteration}.dat"
                cmd = f"cp {join(self.dir_csm, tmp)} {join(self.dir_csm, tmp2)}"
                os.system(cmd)

        # Run Abaqus and check for (licensing) errors
        bool_completed = 0
        attempt = 0
        while not bool_completed and attempt < 10000:
            attempt += 1
            if attempt > 1:
                tools.print_info(f"Warning attempt {attempt-1} in AbaqusSolver failed, new attempt in one minute",
                                 layout='warning')
                time.sleep(60)
                tools.print_info(f"Starting attempt {attempt}")
            if self.timestep == 1:
                cmd1 = f"export PBS_NODEFILE=AbaqusHosts.txt && unset SLURM_GTIDS"
                cmd2 = f"abaqus job=CSM_Time{self.timestep} input=CSM_Time{self.timestep - 1}" \
                    f" cpus={self.cores} user=usr.f output_precision=full interactive >> AbaqusSolver.log 2>&1"
                commands = [cmd1, cmd2]
                self.run_shell(self.dir_csm, commands, name='Abaqus_Calculate')
            else:
                cmd1 = f"export PBS_NODEFILE=AbaqusHosts.txt && unset SLURM_GTIDS"
                cmd2 = f"abaqus job=CSM_Time{self.timestep} oldjob=CSM_Time{self.timestep - 1} input=CSM_Restart" \
                    f" cpus={self.cores} user=usr.f output_precision=full interactive >> AbaqusSolver.log 2>&1"
                commands = [cmd1, cmd2]
                self.run_shell(self.dir_csm, commands, name='Abaqus_Calculate')

            # Check log for completion and or errors
            cmd = "tail -n 10 AbaqusSolver.log > Temp_log.coco"
            self.run_shell(self.dir_csm, [cmd], name='Temp_log')
            log_tmp = os.path.join(self.dir_csm, "Temp_log.coco")
            bool_lic = 1
            with open(log_tmp, "r") as fp:
                for line in fp:
                    if any(x in line for x in ["Licensing error", "license error", "Error checking out Abaqus license"]):
                        bool_lic = 0
            if not bool_lic:
                tools.print_info("Abaqus licensing error", layout='fail')
            elif "COMPLETED" in line:  # Check final line for completed
                bool_completed = 1
            elif bool_lic:  # Final line did not contain "COMPLETED" but also no licensing error detected
                raise RuntimeError("Abaqus did not COMPLETE, unclassified error, see AbaqusSolver.log for extra information")

            # Append additional information to log file
            cmd = f"tail -n 23 CSM_Time{self.timestep}.msg | head -n 15 | sed -e \'s/^[ \\t]*//\' >> AbaqusSolver.log 2>&1"
            self.run_shell(self.dir_csm, [cmd], name='Append_log')

        # Write Abaqus output
        cmd = f"abaqus ./GetOutput.exe CSM_Time{self.timestep} 1 >> AbaqusSolver.log 2>&1"
        self.run_shell(self.dir_csm, [cmd], name='GetOutput')

        # Read Abaqus output data
        for dct in self.interface_output.parameters:
            mp_name = dct['model_part']

            # read in output file for surface nodes
            mp_id = self.model_part_surface_ids[mp_name]
            tmp = f"CSM_Time{self.timestep}Surface{mp_id}Output.dat"
            file_name = join(self.dir_csm, tmp)
            data = np.loadtxt(file_name, skiprows=1)

            # copy output data for debugging
            if self.debug:
                tmp2 = f"CSM_Time{self.timestep}Surface{mp_id}Output_Iter{self.iteration}.dat"
                cmd = f"cp {file_name} {join(self.dir_csm,tmp2)}"
                os.system(cmd)

            if data.shape[1] != self.dimensions:
                raise ValueError(f"given dimension does not match coordinates")

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
        if self.timestep and (self.timestep-1) % self.settings['save_iterations']:
            to_be_removed_suffix = [".com", ".dat", ".mdl", ".msg", ".odb", ".prt", ".res", ".sim", ".sta", ".stt",
                                    "Surface0Cpu0Input.dat", "Surface0Output.dat"]
            cmd = []
            for suffix in to_be_removed_suffix:
                cmd.append(f"rm CSM_Time{self.timestep - 1}{suffix}")
            self.run_shell(self.dir_csm, cmd, name="Remove_previous")

    def finalize(self):
        super().finalize()

    def get_interface_input(self):
        return self.interface_input

    def set_interface_input(self):
        Exception("This solver interface provides no mapping.")  # TODO: remove?

    def get_interface_output(self):
        return self.interface_output

    def set_interface_output(self):
        Exception("This solver interface provides no mapping.")  # TODO: remove?

    def make_elements(self, face_file, output_file):
        firstLoop = 1
        # element = 0
        element_prev = -1
        # point = 0
        point_prev = -1
        element_0 = -1
        point_0 = -1
        count = 0
        element_str = ""

        with open(face_file, 'r') as file:
            for line in file:
                values = line.strip().split()
                element = int(values[0])
                point = int(values[1])
                if element == element_0 and point == point_0:
                    break
                if element == element_prev:
                    if point == point_prev + 1:
                        point_prev = point
                    else:
                        raise ValueError(f"loadpoint number increases by more than one per line for element {element}")
                else:
                    if point == 1:
                        point_prev = point
                        element_prev = element
                        element_str += str(element) + "\n"
                        count += 1
                        if firstLoop:  # Faces contains all values multiple times, but we only want it once
                            element_0 = element
                            point_0 = point
                            firstLoop = 0
                    else:
                        raise ValueError(f"Load point number does not start at 1 for element {element}")

        element_str = f"{count}\n{point_prev}\n" + element_str
        with open(output_file, "w") as file:
            file.write(element_str)

    def run_shell(self, work_dir, commands, wait=True, name='script', delete=True):
        script = f'{work_dir}/{name}.sh'
        with open(script, 'w') as file:
            file.write('#!/bin/bash\n')
            file.write(f'cd {work_dir}\n')
            for line in commands:
                file.write(line + '\n')
        os.chmod(script, 0o700)
        if wait:
            p = subprocess.call(script, shell=True)
        else:
            p = subprocess.Popen(script, shell=True)
        if delete:
            os.system("rm " + script)
        return p

    def replace_fortran(self, line, orig, new):
        """The length of a line in FORTRAN 77 is limited, replacing working directories can exceed this limit.
        This functions splits these strings over multiple lines."""

        ampersand_location = 6
        char_limit = 72

        if "|" in line:
            temp = line.replace(orig, new)
            N = len(temp)

            if N > char_limit:
                count = 0
                line = ""
                line += temp[0:char_limit] + "\n"
                count += char_limit
                while count < N:
                    temp_string = temp[count:count + char_limit - ampersand_location]
                    n = len(temp_string)
                    count += n
                    if count < N:  # need to append an additional new line
                        line += "     &" + temp_string + "\n"
                    else:
                        line += "     &" + temp_string
            else:
                line = temp

        return line

    def write_start_and_restart_inp(self, input_file, output_file, restart_file):
        bool_restart = 0

        rf = open(restart_file, "w")
        of = open(output_file, "w")

        rf.write("*HEADING \n")
        rf.write("*RESTART, READ \n")

        with open(input_file) as f:
            line = f.readline()
            while line:
                if "*step" in line.lower():
                    contents = line.split(",")  # Split string on commas
                    line_new = ''
                    for s in contents:
                        if s.strip().startswith("inc="):
                            numbers = re.findall("\d+", s)
                            if (not self.subcycling) and int(numbers[0]) != 1:
                                raise NotImplementedError(f"inc={numbers[0]}: subcycling was not requested but maxNumInc > 1.")
                            else:
                                line_new += f" inc={self.maxNumInc},"
                        else:
                            line_new += s+","
                    line_new = line_new[:-1]+"\n"  # Remove the last added comma and add a newline

                    of.write(line_new)
                    if bool_restart:
                        rf.write(line_new)
                    line = f.readline()
                elif "*dynamic" in line.lower():
                    contents = line.split(",")
                    for s in contents:
                        if "application" in s.lower():
                            contents_B = s.split("=")
                            app = contents_B[1].lower().strip()
                            if app == "quasi-static" or app == "moderate dissipation":
                                if not self.subcycling:
                                    line_2 = f"{self.delta_t}, {self.delta_t},\n"
                                else:
                                    line_2 = f"{self.initialInc}, {self.delta_t}, {self.minInc}, {self.maxInc}\n"
                            else:
                                raise NotImplementedError(
                                    f"{contents_B[1]} not available: Currently only quasi-static and moderate dissipation are implemented for the Abaqus wrapper")
                    of.write(line)
                    if bool_restart:
                        rf.write(line)
                        f.readline()  # need to skip the next line
                        of.write(line_2)  # Change the time step in the Abaqus step
                        if bool_restart:
                            rf.write(line_2)  # Change the time step in the Abaqus step (restart-file)
                        line = f.readline()
                    elif "*static" in line.lower():
                        of.write(line)
                        if bool_restart:
                            rf.write(line)
                        f.readline()  # need to skip the next line
                        if not self.subcycling:
                            raise NotImplementedError(
                                "Only Static with subcycling is implemented for the Abaqus wrapper")
                        line_2 = f"{self.initialInc}, {self.delta_t}, {self.minInc}, {self.maxInc}\n"
                        f.readline()
                    of.write(line_2)  # Change the time step in the Abaqus step
                    if bool_restart:
                        rf.write(line_2)  # Change the time step in the Abaqus step (restart-file)
                    line = f.readline()
                else:
                    of.write(line)
                    if bool_restart:
                        rf.write(line)
                    line = f.readline()
                if "** --"in line:
                    bool_restart = 1
        rf.close()
        of.close()

    def write_loads(self):
        for dct in self.interface_input.parameters:
            mp_name = dct['model_part']
            mp_id = self.model_part_surface_ids[mp_name]
            model_part = self.model.get_model_part(mp_name)

            pressure = self.interface_input.get_variable_data(mp_name, 'pressure')
            traction = self.interface_input.get_variable_data(mp_name, 'traction')
            data = np.hstack((pressure, traction[:, :self.dimensions]))
            fmt = (self.dimensions + 1) * '%27.17e'  # format of load input file should be exactly this for USR.f
            tmp = f'CSM_Time{self.timestep}Surface{mp_id}Cpu0Input.dat'
            file_name = join(self.dir_csm, tmp)
            np.savetxt(file_name, data, fmt=fmt, header=f'{model_part.size}', comments='')

            # Start of a simulation with ramp, needs an initial load at time 0: set at zero load
            if self.iteration == 1 and self.timestep == 1 and self.ramp:
                tmp = f'CSM_Time{self.timestep-1}Surface{mp_id}Cpu0Input.dat'
                file_name = join(self.dir_csm, tmp)
                np.savetxt(file_name, data * 0.0, fmt=fmt, header=f'{model_part.size}', comments='')
