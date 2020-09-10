import numpy as np
import os
import linecache
import sys
import time
import copy
import subprocess

from coconut import data_structure
from coconut.coupling_components.component import Component
from coconut.coupling_components.interface import Interface


def Create(parameters):
    return SolverWrapperOpenFOAM_41(parameters)


class SolverWrapperOpenFOAM_41(Component):
    def __init__(self, parameters):
        super().__init__()
        
        #Settings
        self.settings = parameters["settings"]
        self.working_directory = self.settings["working_directory"].GetString()
#         input_file = self.settings["input_file"].GetString()
#         settings_file_name = os.path.join(working_directory, input_file)
#         with open(settings_file_name, 'r') as settings_file:
#             self.settings.AddParameters(cs_data_structure.Parameters(settings_file.read()))
        self.moduleType="OpenFOAM" #hard-coded, cannot be altered by naive user
        self.moduleVersion="4.1" #hard-coded, cannot be altered by naive user
        self.module=self.moduleType+"/"+self.moduleVersion
        self.application=self.settings["application"].GetString() #What type of OF-solver to be used - solver requires adaptation before running with OpenFOAM - name of adapted solver starts with 'CoCoNuT_'
        self.dimensions=self.settings["dimensions"].GetInt()
        if (self.dimensions != 2) and (self.dimensions != 3):
            sys.exit("OpenFOAM-case should be 2D or 3D.")
        self.dt = self.settings["dt"].GetDouble()  # Time step size
        self.start_time=self.settings["start_time"].GetDouble() # Start time - also the name of folder containing data at this time
        self.end_time=self.settings["end_time"].GetDouble() # End time
        self.cores=self.settings["cores"].GetInt() # Number of cores to be used in the OpenFOAM-calculation
        self.decomposeMethod=self.settings["decomposeMethod"].GetString() #Decomposition-method, can be "simple", "scotch" 
        self.newtonmax = self.settings["newtonmax"].GetInt()  # Maximal number of Newton iterations
        self.newtontol = self.settings["newtontol"].GetDouble()  # Tolerance of Newton iterations
        self.write_interval = self.settings["write_interval"].GetInt() # Number of time steps between consecutive saves performed by OpenFOAM 
        self.write_precision = self.settings["write_precision"].GetInt() # writePrecision-parameter in OpenFOAM
        self.time_precision = self.settings["time_precision"].GetInt() # timePrecision-parameter in OpenFOAM
        self.time_format = "fixed"
        self.boundary_names = [_.GetString() for _ in self.settings['boundary_names'].list()] # boundary_names is the set of boundaries where the moving interface is located (will be used to go through OF-files)
        self.meshmotion_solver=self.settings["meshmotion_solver"].GetString()
        self.diffusivity=self.settings["diffusivity"].GetString()
        
        #Check that the boundary_names and the interface_input and interface_output are defined consistently in the JSON-file
        #For every boundary_name element, there should be one interface_input (boundary_name+"_input") element and one interface_output (boundary_name+"_output") element.
        #Make the distinction between both: interface_input/output are names of the pyKratos ModelPart - boundary_names is the name of the boundary as defined in OpenFOAM!
        if len(self.boundary_names) != len(self.settings['interface_input'].keys()):
            sys.exit("Interface_input and boundary_names should have the same length; one interface_input element corresponds with one boundary_name element!")
        if len(self.boundary_names) != len(self.settings['interface_output'].keys()):
            sys.exit("Interface_output and boundary_names should have the same length; one interface_output element corresponds with one boundary_name element!")
        index=0
        for boundary in self.boundary_names:
            key_input = self.settings['interface_input'].keys()[index]
            if not("_input" in key_input):
                sys.exit('Please make that all interface_input elements correspond to a boundary_names element followed by "_input".')
            else:
                keyBoundary=key_input.replace("_input","")
                if keyBoundary != boundary:
                    sys.exit('For OpenFOAM, please make sure that every boundary_names element is linked to an interface_input element with the following name: <boundary_names element>"_input". The corresponding elements in boundary_names and interface_input should have the same index in their respective arrays!')
            key_output = self.settings['interface_output'].keys()[index]
            if not("_output" in key_output):
                sys.exit('Please make that all interface_output elements correspond to a boundary_names element followed by "_output".')
            else:
                keyBoundary=key_output.replace("_output","")
                if keyBoundary != boundary:
                    sys.exit('For OpenFOAM, please make sure that every boundary_names element is linked to an interface_output element with the following name: <boundary_names element>"_output". The corresponding elements in boundary_names and interface_output should have the same index in their respective arrays!')
                if key_output.replace("_output","_input") != key_input:
                    sys.exit("Please make sure that the interface_input and interface_output elements occur in such a way that the corresponding elements have the same index.")
            index+=1
 
        #Check that the correct modules have been loaded
        self.check_software()  
        
        #Remove possible CoCoNuT-message from previous interrupt
        self.remove_all_messages()
        
        # Creating OpenFOAM-files - raw dictionary files are predefined in the solver_wrapper folder (and should not be moved)
        # DecomposeParDict: replace raw settings by actual settings defined by user in json-file
        if self.cores > 1: #Only if calculating in parallel
            decomposeParDict_raw_name=os.path.join(os.path.realpath(os.path.dirname(__file__)),"decomposeParDict_raw")
            decomposeParDict_name=os.path.join(self.working_directory,"system/decomposeParDict")
            with open(decomposeParDict_raw_name,'r') as rawFile:
                with open(decomposeParDict_name,'w') as newFile:
                    for line in rawFile:
                        line=line.replace('|CORES|',str(self.cores))
                        line=line.replace('|DECOMPOSEMETHOD|',str(self.decomposeMethod))
                        newFile.write(line)
            rawFile.close()
            newFile.close()
            self.write_footer(decomposeParDict_name) 
            # # OpenFOAM-fields are decomposed automatically if you work in parallel
            # os.system("cd " + self.working_directory + "; decomposePar -force -time "+ str(self.start_time) + " &> log.decomposePar;")

        # ControlDict: replace raw settings by actual settings defined by user in json-file AND add function objects to write pressure and wall shear stress
        controlDict_raw_name=os.path.join(os.path.realpath(os.path.dirname(__file__)),"controlDict_raw")
        controlDict_name=os.path.join(self.working_directory,"system/controlDict")
        with open(controlDict_raw_name,'r') as rawFile:
            with open(controlDict_name,'w') as newFile:
                for line in rawFile:
                    line=line.replace('|APPLICATION|',str(self.application))
                    line=line.replace('|START_TIME|',str(self.start_time))
                    line=line.replace('|END_TIME|',str(self.end_time))
                    line=line.replace('|DT|',str(self.dt))
                    line=line.replace('|WRITE_INTERVAL|',str(self.write_interval))
                    line=line.replace('|WRITE_PRECISION|',str(self.write_precision))
                    line=line.replace('|TIME_PRECISION|',str(self.time_precision))
                    if '|BOUNDARY_NAMES|' in line:
                        firstBoundary=True
                        for interfaces in self.boundary_names:
                            if firstBoundary:
                                boundary_name_temp = "(" + interfaces
                                firstBoundary=False
                            else:
                                boundary_name_temp += " , " + interfaces
                        boundary_name_temp += ")"                          
                        line=line.replace('|BOUNDARY_NAMES|',boundary_name_temp)
                    newFile.write(line)
        rawFile.close()
        newFile.close()
        nKey=0
        if len(self.boundary_names) == 1:
            for key in self.boundary_names:
                self.write_controlDict_function(controlDict_name,"surfaceRegion","libfieldFunctionObjects","PRESSURE",key,True,False)
                self.write_controlDict_function(controlDict_name,"wallShearStress","libfieldFunctionObjects","wallShearStress",key,False,False)
                self.write_controlDict_function(controlDict_name,"surfaceRegion","libfieldFunctionObjects","TRACTION",key,False,True)
        else:
            for key in self.boundary_names:
                if nKey == 0:
                    self.write_controlDict_function(controlDict_name,"surfaceRegion","libfieldFunctionObjects","PRESSURE",key,True,False)
                    self.write_controlDict_function(controlDict_name,"wallShearStress","libfieldFunctionObjects","wallShearStress",key,False,False)
                    self.write_controlDict_function(controlDict_name,"surfaceRegion","libfieldFunctionObjects","TRACTION",key,False,False)
                else:
                    self.write_controlDict_function(controlDict_name,"surfaceRegion","libfieldFunctionObjects","PRESSURE",key,False,False)
                    if nKey == (len(self.boundary_names)-1):
                        self.write_controlDict_function(controlDict_name,"surfaceRegion","libfieldFunctionObjects","TRACTION",key,False,True)
                    else:
                        self.write_controlDict_function(controlDict_name,"surfaceRegion","libfieldFunctionObjects","TRACTION",key,False,False)
                nKey += 1
        self.write_footer(controlDict_name)
        # DynamicMeshDict: replace raw settings by actual settings defined by user in json-file 
        dynamicMeshDict_raw_name=os.path.join(os.path.realpath(os.path.dirname(__file__)),"dynamicMeshDict_raw")
        dynamicMeshDict_name=os.path.join(self.working_directory,"constant/dynamicMeshDict")
        strBoundary=""
        for boundary in self.boundary_names:
            strBoundary=strBoundary+" "+str(boundary)
        with open(dynamicMeshDict_raw_name,'r') as rawFile:
            with open(dynamicMeshDict_name,'w') as newFile:
                for line in rawFile:
                    line=line.replace('|MESHMOTION_SOLVER|',str(self.meshmotion_solver))
                    line=line.replace('|DIFFUSIVITY|',str(self.diffusivity))
                    line=line.replace('|NUM_INTERFACE_INPUT|',str(len(self.settings['boundary_names'].list())))
                    line=line.replace('|INTERFACE_INPUT|',strBoundary)
                    newFile.write(line)
        rawFile.close()
        newFile.close()
        nKey=0
        self.write_footer(dynamicMeshDict_name)
        
        # Creating Model
        self.model = data_structure.Model()
        print("The model for OpenFOAM will be created. Please make sure all patch names given under the 'interface' setting are also found in the mesh used in OpenFOAM (see 'constant/polyMesh') \n")
        
        # Creating ModelParts and adding variables to these ModelParts - should happen before node addition
        for key,value in (self.settings['interface_input'].items()+self.settings['interface_output'].items()):
            self.model.CreateModelPart(key)
            mp=self.model[key]
            for var_name in value.list():
                var=vars(data_structure)[var_name.GetString()]
                mp.AddNodalSolutionStepVariable(var)
            
        # Adding nodes to ModelParts - should happen after variable definition; writeCellcentres writes cellcentres in internal field and face centres in boundaryField
        os.system("cd "+ self.working_directory + "; writeCellCentres -time " + str(self.start_time) + " &> log.writeCellCentres;")
        nKey=0
        for boundary in self.boundary_names:
            source_file = self.working_directory + "/constant/polyMesh"
            node_ids, node_coords, face_centres = self.Get_Point_IDs(boundary,source_file)

            mp_input = self.model[boundary + "_input"]
            for i in np.arange(0,len(node_ids)):
                mp_input.CreateNewNode(node_ids[i],node_coords[i, 0], node_coords[i, 1], node_coords[i, 2])

            index = 0
            mp_output = self.model[boundary + "_output"]
            for i in np.arange(0, len(face_centres)):
                mp_output.CreateNewNode(index, face_centres[i, 0], face_centres[i, 1], face_centres[i, 2])
                index+=1


        # Create CoSimulationInterfaces
        self.interface_input = Interface(self.model, self.settings["interface_input"])
        self.interface_output = Interface(self.model, self.settings["interface_output"])

        # Create Variables
        self.pressure=vars(data_structure)['PRESSURE']
        self.shear=vars(data_structure)['TRACTION']
        self.displacement=vars(data_structure)['DISPLACEMENT']
        
             
    def Initialize(self):
        super().Initialize()
                
        # Define timestep and physical time
        self.timestep=0
        self.physical_time=self.start_time
        
        # If no pointDisplacement file is defined yet, initialize a pointDisplacement file in the start time folder
        # Normally, after restart or from the second iteration onwards, a pointDisplacement-file already exists. In that case, that pointDisplacement-file will be used (and is NOT overwritten)
        pointDisp_raw_name=os.path.join(os.path.realpath(os.path.dirname(__file__)),"pointDisplacement_raw")
        pointDisp_name=os.path.join(self.working_directory,str(self.physical_time),"pointDisplacement")
        if not(os.path.isfile(pointDisp_name)):
            self.write_pointDisplacement_file(pointDisp_raw_name,pointDisp_name)
        # if self.cores == 1:
        #     pointDisp_name=os.path.join(self.working_directory,str(self.physical_time),"pointDisplacement")
        #     if not(os.path.isfile(pointDisp_name)):
        #         self.write_pointDisplacement_file(pointDisp_raw_name,pointDisp_name)
        # else:
        #     for p in np.arange(self.cores):
        #         pointDisp_name=os.path.join(self.working_directory,"processor"+str(p),str(self.physical_time),"pointDisplacement")
        #         if not(os.path.isfile(pointDisp_name)):
        #             self.write_pointDisplacement_file(pointDisp_raw_name,pointDisp_name)

        ## Do a first decomposition
        if self.cores > 1:
            os.system("cd " + self.working_directory + "; decomposePar -force -time "+ str(self.start_time) + " &> log.decomposePar;")
        
        # Don't forget to start OpenFOAM-loop!
        if self.cores == 1:
            cmd = self.application + "&> log." + self.application
        else:
            cmd = "mpirun -np " + str(self.cores) + " " + self.application + " -parallel &> log." + self.application

        self.openfoam_process = subprocess.Popen(cmd, cwd=self.working_directory, shell=True) 

        ### CoConuT_OpenFOAMSolver is running from here on
                                        

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()

        # For parallel runs need to create a folder with the correct time stamp for decomposition of pointDisplacement_Next
        # For serial runs, this folder will normally be present
        timestamp = '{:.{}f}'.format(self.physical_time, self.time_precision)
        path = os.path.join(self.working_directory, timestamp)
        if self.cores > 1:
            os.system('mkdir '+path)
        elif self.physical_time == 0: # for serial also need to make a folder 0.0000 with specified precision
            os.system('mkdir '+path)

        # # The displacement of the FSI interface is passed through pointDisplacement_Next, which is prepared here
        # pointDisp_raw_name=os.path.join(os.path.realpath(os.path.dirname(__file__)),"pointDisplacement_raw")
        # pointDisp_name=os.path.join(self.working_directory,str(self.physical_time),timestamp,"pointDisplacement_Next")
        # self.write_pointDisplacement_file(pointDisp_raw_name,pointDisp_name)

        # Prepare new time step folder and reset the number of iterations
        self.timestep += 1
        self.iteration = 0
        self.physical_time += self.dt

        self.prev_timestamp = timestamp
        self.cur_timestamp = '{:.{}f}'.format(self.physical_time, self.time_precision)

        if self.cores <2: #if serial
            newPath=os.path.join(self.working_directory, self.cur_timestamp)
            if os.path.isdir(newPath):
                print("\n\n\n Warning! In 5s, CoCoNuT will overwrite existing time step folder: "+str(newPath) + ". \n\n\n")
                time.sleep(5)
                os.system("rm -rf " + newPath)
            os.system("mkdir "+ newPath)
        else :
            for i in np.arange(self.cores):
                newPath=os.path.join(self.working_directory, "processor"+str(i), self.cur_timestamp)
                if os.path.isdir(newPath):
                    if i==0:
                        print("\n\n\n Warning! In 5s, CoCoNuT will overwrite existing time step folder in processor-subfolders. \n\n\n")
                        time.sleep(5)
                    os.system("rm -rf " + newPath)
                os.system("mkdir "+ newPath)
        print('\t Time step '+str(self.timestep))
        
        self.send_message('next') # Let OpenFOAM go to next time step
        self.wait_message('next_ready') # Let OpenFOAM wait for input data
    

    def SolveSolutionStep(self, interface_input):
        self.iteration += 1
        print(f'\t\tIteration {self.iteration}')

        # store incoming displacements
        self.interface_input.SetPythonList(interface_input.GetPythonList())

        # update X,Y,Z in interface
        for key in [_[0] for _ in self.interface_input.model_parts_variables]:
            for node in self.model[key].Nodes:
                disp = node.GetSolutionStepValue(self.displacement)
                node.X = node.X0 + disp[0]
                node.Y = node.Y0 + disp[1]
                node.Z = node.Z0 + disp[2]

        # write interface data to OpenFOAM-file
        self.write_node_input()
            
        # let Fluent run, wait for data
        self.send_message('continue')
        self.wait_message('continue_ready')

        # read data from OpenFOAM
        self.read_node_output()
        
        # return interface_output object
        return self.interface_output.deepcopy()


    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        # Let OpenFOAM check whether it needs to save this timestep (in OF-solver: runTime.write())
        
        if not(self.timestep % self.write_interval):
            self.send_message('save')
            self.wait_message('save_ready')
#         else:
#             # This is done to remove the OpenFOAM-subfolder containing the wallShearStress, pressure and displacement for that particular time step
#             os.system("rm -r " + os.path.join(self.working_directory, self.physical_time))
#             pass     
            
            
    def Finalize(self):
        super().Finalize()
        
        self.send_message('stop')
        self.wait_message('stop_ready')
        
        self.openfoam_process.kill()
                
        print("OpenFOAM was stopped with the Finalize() method defined in CoCoNuT.")


    def GetInterfaceInput(self):
        return self.interface_input.deepcopy()


    def SetInterfaceInput(self):
        Exception("This solver interface provides no mapping.")


    def GetInterfaceOutput(self):
        return self.interface_output.deepcopy()


    def SetInterfaceOutput(self):
        Exception("This solver interface provides no mapping.")

    
    def write_header(self,fileLoc,className,objectName):
        f=open(fileLoc,'w')
        f.write(r'/*--------------------------------*- C++ -*----------------------------------*\\'+"\n")
        f.write(r'| =========                 |                                                 |'+"\n")
        f.write(r'| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |'+"\n")
        f.write(r'|  \\    /   O peration     | Version:  4.x                                   |'+"\n")
        f.write(r'|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |'+"\n")
        f.write(r'|    \\/     M anipulation  |                                                 |'+"\n")
        f.write(r'\*---------------------------------------------------------------------------*/'+"\n")
        f.write(r'FoamFile'+"\n")
        f.write(r'{'+"\n")
        f.write('\t version \t\t 4.1;'+"\n")
        f.write('\t format \t\t ascii;'+"\n")
        f.write('\t class \t\t ' + className + ';'+"\n")
        f.write('\t object \t\t ' + objectName + ';'+"\n")
        f.write('}'+"\n")
        f.write(r'// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //'+"\n")
        f.write("\n")  
        f.close()
        
    def read_node_output(self):

        nKey=0
        for boundary in self.boundary_names:
            # specify location of pressure and traction
            tractionName="TRACTION_"+boundary
            pressureName="PRESSURE_"+boundary
            mp = self.model[boundary+"_output"]
            nFaces_tot = mp.NumberOfNodes()
            wss_tmp=np.zeros([nFaces_tot,3])
            pres_tmp=np.zeros([nFaces_tot,1])
            wss_file= os.path.join(self.working_directory, "postProcessing", tractionName, "surface", self.cur_timestamp,"wallShearStress_patch_"+boundary+".raw")
            pres_file= os.path.join(self.working_directory, "postProcessing", pressureName, "surface", self.cur_timestamp,"p_patch_"+boundary+".raw")

            # read traction
            counter = 0
            nlines = 0
            while (nlines<nFaces_tot+2) and counter < 100:
                if os.path.isfile(wss_file):
                    nlines = sum(1 for line in open(wss_file))

                time.sleep(0.01)
                counter +=1

            if counter == 100:
                raise RuntimeError("Timed out waiting for wss file: "+wss_file)

            f=open(wss_file,'r')
            fLines=f.readlines()
            index_start=2
            for i in np.arange(nFaces_tot):
                wss_tmp[i,0]=fLines[index_start+i].split()[3]
                wss_tmp[i,1]=fLines[index_start+i].split()[4]
                wss_tmp[i,2]=fLines[index_start+i].split()[5].split("\n")[0]
            f.close()
            # read pressure
            f=open(pres_file,'r')
            fLines=f.readlines()
            index_start=2
            it=0
            for i in np.arange(nFaces_tot):
                val=fLines[index_start+i].split()[3].split("\n")[0]
                pres_tmp[i,0]=float(val)
            f.close()
            # store pressure and traction in Nodes
            index=0
            for node in mp.Nodes:
                node.SetSolutionStepValue(self.shear, 0, wss_tmp[index])
                node.SetSolutionStepValue(self.pressure, 0, pres_tmp[index])
                index += 1
            
            # go to next interface
            nKey += 1
      
        
    #writeFooter: to write OpenFOAM-footer in file at location 'fileLoc'
    def write_footer(self,fileLoc):
        f=open(fileLoc,'a+')
        f.write("\n")
        f.write(r'// ************************************************************************* //'+"\n")
        f.close()        
    
    
    def write_node_input(self):
        # The displacement of the FSI interface is passed through the file "pointDisplacement_Next"
        # This function will prepare that file in a "serial format" and then decompose it for parallel operation

        # (An alternative is to keep track of the mapping between processors, but this entails a lot of "exceptions" in terms of reading the files)

        #self.prev_timestamp
        #self.cur_timestamp

        pointDisp_raw_name = os.path.join(os.path.realpath(os.path.dirname(__file__)), "pointDisplacement_raw")
        pointDisp_name = os.path.join(self.working_directory, self.prev_timestamp, "pointDisplacement_Next")
        self.write_pointDisplacement_file(pointDisp_raw_name, pointDisp_name)

        disp_file = pointDisp_name

        nKey = 0
        for boundary in self.boundary_names:
            mp = self.model[boundary + "_input"]

            startNr = self.find_string_in_file(boundary, disp_file)
            os.system("head -n " + str(startNr + 1) + " " + disp_file + " > tempDisp")

            with open('tempDisp', 'a+') as file:
                file.write("\t { \n")
                file.write("\t\t type  \t fixedValue; \n")
                file.write('\t\t value \t nonuniform List<vector> ( \n')
                for node in mp.Nodes:
                    dispX = node.X - node.X0
                    dispY = node.Y - node.Y0
                    dispZ = node.Z - node.Z0
                    file.write(' (' + f'{dispX:27.17e} {dispY:27.17e} {dispZ:27.17e}' + ') \n')
                file.write(');\n')
            # if self.nNodes_tot < 11:  # 10 or less elements on interface
            #     with open('tempDisp', 'a+') as file:
            #         file.write("\t { \n")
            #         file.write("\t\t type  \t fixedValue; \n")
            #         file.write('\t\t value \t nonuniform List<vector> (')
            #         for node in mp.Nodes:
            #             dispX = node.X - node.X0
            #             dispY = node.Y - node.Y0
            #             dispZ = node.Z - node.Z0
            #             file.write(' (' + f'{dispX:27.17e} {dispY:27.17e} {dispZ:27.17e}' + ') ')
            #         file.write(');\n')
            #     file.close()
            # else:
            #     with open('tempDisp', 'a+') as file:
            #         file.write("\t { \n")
            #         file.write("\t\t type  \t fixedValue; \n")
            #         file.write('\t\t value \t nonuniform List<vector> ( \n')
            #         for node in mp.Nodes:
            #             dispX = node.X - node.X0
            #             dispY = node.Y - node.Y0
            #             dispZ = node.Z - node.Z0
            #             file.write(' (' + f'{dispX:27.17e} {dispY:27.17e} {dispZ:27.17e}' + ') \n')
            #         file.write(');\n')
            #     file.close()
            os.system("wc -l " + disp_file + " > lengthDisp")
            lengthDisp_file = open("lengthDisp", 'r')
            length_disp = int(lengthDisp_file.readline().split(" ")[0])
            lengthDisp_file.close()
            os.system("tail -n " + str(length_disp - (startNr + 1)) + " " + disp_file + " > tempDisp2")
            startToEndNr = self.find_string_in_file("}", "tempDisp2")
            os.system("tail -n " + str(length_disp - (startNr + 1) - startToEndNr) + " " + disp_file + " > tempDisp3")
            os.system("cat tempDisp tempDisp3 > " + disp_file)
            os.system("rm tempDisp* lengthDisp")
            nKey += 1


        #     ### Timestep zero should be reconstructPar with the option -withZero"
        # os.system(f"reconstructPar -time '{self.physical_time-self.dt:.6f}' -fields '(pointDisplacement_Next)'")
        # nKey=0
        # for boundary in self.boundary_names:
        #     mp = self.model[boundary+"_input"]
        #     if self.cores<2: #serial run
        #         disp_file = os.path.join(self.working_directory, f"{self.physical_time-self.dt:.6f}", 'pointDisplacement_Next') #Need to update pointdisplacement in the previous time directory
        #     else:
        #         for c in range(0,self.cores):
        #             disp_file = os.path.join(self.working_directory,f"processor{c}", f"{self.physical_time - self.dt:.6f}",'pointDisplacement_Next')
        #
        #     startNr=self.find_string_in_file(boundary, disp_file)
        #     os.system("head -n " + str(startNr+1) + " " + disp_file + " > tempDisp")
        #     if self.nNodes_tot < 11: # 10 or less elements on interface
        #         with open('tempDisp', 'a+') as file:
        #             file.write("\t { \n")
        #             file.write("\t\t type  \t fixedValue; \n")
        #             file.write('\t\t value \t nonuniform List<vector> (')
        #             for node in mp.Nodes:
        #                 dispX=node.X-node.X0
        #                 dispY=node.Y-node.Y0
        #                 dispZ=node.Z-node.Z0
        #                 file.write(' (' + f'{dispX:27.17e} {dispY:27.17e} {dispZ:27.17e}'+') ')
        #             file.write(');\n')
        #         file.close()
        #     else :
        #         with open('tempDisp', 'a+') as file:
        #             file.write("\t { \n")
        #             file.write("\t\t type  \t fixedValue; \n")
        #             file.write('\t\t value \t nonuniform List<vector> ( \n')
        #             for node in mp.Nodes:
        #                 dispX=node.X-node.X0
        #                 dispY=node.Y-node.Y0
        #                 dispZ=node.Z-node.Z0
        #                 file.write(' (' + f'{dispX:27.17e} {dispY:27.17e} {dispZ:27.17e}'+') \n')
        #             file.write(');\n')
        #         file.close()
        #
        #     os.system("wc -l " + disp_file + " > lengthDisp")
        #     lengthDisp_file=open("lengthDisp",'r')
        #     length_disp=int(lengthDisp_file.readline().split(" ")[0])
        #     lengthDisp_file.close()
        #     os.system("tail -n " + str(length_disp-(startNr+1)) + " " + disp_file + " > tempDisp2")
        #     startToEndNr=self.find_string_in_file("}", "tempDisp2")
        #     os.system("tail -n " + str(length_disp-(startNr+1)-startToEndNr) + " " + disp_file + " > tempDisp3")
        #     os.system("cat tempDisp tempDisp3 > "+ disp_file)
        #     os.system("rm tempDisp* lengthDisp")
        #     nKey += 1

        if self.cores > 1:
            os.system("cd " + self.working_directory + "; decomposePar -fields -time " + self.prev_timestamp + " &> log.decomposePar;")


    def write_controlDict_function(self, filename, funcname, libname, varname, patchname, writeStart, writeEnd):
        with open(filename,'a+') as file:
            if writeStart:
                file.write("functions \n")
                file.write("{ \n ")
            if varname=="wallShearStress":
                file.write(" \n \t " + varname + " \n")
            else:
                file.write(" \n \t " + varname + "_" + patchname +" \n")
            file.write("\t { \n")
            file.write("\t\t type  \t " + funcname + "; \n")
            file.write('\t\t libs \t ("' + libname + '.so"); \n')
            file.write('\t\t executeControl \t timeStep; \n')
            file.write('\t\t executeInterval \t 1; \n')
            file.write('\t\t writeControl \t timeStep; \n')
            file.write('\t\t writeInterval \t 1; \n')
            file.write('\t\t timeFormat \t fixed; \n')
            file.write(f'\t\t timePrecision \t {self.time_precision}; \n')
            if funcname=="surfaceRegion":
                file.write('\t\t operation \t none; \n')
                file.write('\t\t writeFields \t true; \n')
                file.write('\t\t surfaceFormat \t raw; \n')
                file.write('\t\t regionType \t patch; \n')
                file.write('\t\t name \t ' + patchname + ' ; \n')
                file.write('\t\t fields \n')
                file.write('\t\t ( \n')
                if varname == "PRESSURE":
                    file.write('\t\t\t p \n ')
                elif varname == "TRACTION":
                    file.write('\t\t\t wallShearStress \n')
                file.write("\t\t ); \n")
            elif funcname=="wallShearStress":
#                 file.write('\t\t patches ( ' + patchname + ' ); \n')
                file.write('\t\t log \t false; \n')
            file.write("\t } \n\n")
            if writeEnd:
                file.write("} \n ")
            if varname == "PRESSURE":
                print("\n\n Please check the 'rho' option in the static pressure definition in controlDict! This might vary from OF-solver to OF-solver.\n\n")

        file.close()
            
    
    def write_pointDisplacement_file(self,pointDisp_raw_name,pointDisp_name):
        with open(pointDisp_raw_name,'r') as rawFile: 
                with open(pointDisp_name,'w') as newFile:
                    for line in rawFile:
                        newFile.write(line)
                    nKey=0
                    for boundary in self.boundary_names: 
                        newFile.write(" \n \t "  + boundary +" \n")
                        newFile.write("\t { \n")
                        newFile.write("\t\t type  \t fixedValue; \n")
                        newFile.write("\t\t value \t uniform (0 0 0); \n")
                        newFile.write("\t } \n")
                    newFile.write("} \n")
                    newFile.close()
                    self.write_footer(pointDisp_name)
        rawFile.close()

    
    def find_string_in_file(self,string_name,file_name):
        index=-1
        with open(file_name) as f:
            for num, line in enumerate(f):
                if string_name in line:
                    index=num
                    break
        f.close()
        return index
    
    def send_message(self, message):
        file = os.path.join(self.working_directory, message + ".coco")
        open(file, 'w').close()
        return


    def wait_message(self, message):
        waitTimeLimit=10*60 # 10 minutes maximum waiting time for a single flow solver iteration
        cumulTime=0
        file = os.path.join(self.working_directory, message + ".coco")
        while not os.path.isfile(file):
            time.sleep(0.01)
            cumulTime += 0.01
            if cumulTime > waitTimeLimit:
                os.system("pkill " + self.application)
                sys.exit("CoCoNuT timed out in the OpenFOAM solver_wrapper, waiting for message: "+ message + ".coco.")
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
        if os.system(self.application+' -help &> checkSoftware') != 0:
            sys.exit("You either did not load the module for OpenFOAM/4.1, did not compile the solver you are trying to use or did not source $FOAM_BASH prior to execution of CoCoNuT.")

        # The statement above could work for other version of OpenFOAM is the solver is also compiled for that version. Therefore, the version number is checked explicity (don't forget to remove the end-of-line variable at the end of the String versionNr
        with open('checkSoftware','r') as f:
            lastLine=f.readlines()[-2] # Second last line contains 'Build: XX' with XX the version number 
        f.close()
        os.system('rm checkSoftware')
        versionNr=lastLine.split(' ')[-1]
        if versionNr[:-1] != self.moduleVersion :
                sys.exit("OpenFOAM 4.1 should be loaded! Currently, another version of OpenFOAM is loaded")

    # def read_boundary_file(self,source_file,boundary):
    #     try:
    #         lineNameNr = self.find_string_in_file(boundary, source_file)
    #         lineStartNr = lineNameNr + 7  # In case of non-uniform list, this is where the list of values in the sourceFile starts
    #         nNodesIndex = lineNameNr + 5  # On this line, the number of elements in the boundary
    #         os.system("awk NR==" + str(nNodesIndex) + " " + source_file + " > nNodes")
    #         nNodes_file = open("nNodes", 'r')
    #         nNodes = int(nNodes_file.readline())
    #         nNodes_file.close()
    #         os.system("rm nNodes")
    #         tempCoordFile = np.ones([nNodes, 1]) * float("inf")
    #         for j in np.arange(nNodes):
    #             tempCoordFile[j, 0] = float(linecache.getline(source_file, lineStartNr + j))
    #     except ValueError:  # If ValueError is triggered, it means that the source-file has a uniform coordinate in the axis you are currently looking
    #         if not 'nNodes' in locals():  # If the first coordinate-file has a uniform value, the variable 'nNodes'
    #             # does not exist, so you should check whether this variable exists
    #             check_file = self.working_directory + "/" + str(
    #                 self.start_time) + "/ccy"  # if 'nNodes' does not exist, read the second sourceFile to know the number of rows
    #             lineNameNr = self.find_string_in_file(boundary, check_file)
    #             nNodesIndex = lineNameNr + 5  # On this line, the number of cell centers on the inlet is stated
    #             os.system("awk NR==" + str(nNodesIndex) + " " + check_file + " > nNodes")
    #             nNodes_file = open("nNodes", 'r')
    #             nNodes = int(nNodes_file.readline())
    #             nNodes_file.close()
    #             os.system("rm nNodes")
    #         indexUV = lineNameNr + 4
    #         os.system("awk NR==" + str(indexUV) + " " + source_file + " > unifValue")
    #         unifValue_file = open("unifValue", 'r')
    #         unifValue = float(unifValue_file.readline().split()[-1][
    #                           0:-1])  # First '-1' makes sure the value is read, but this still contains a
    #         # semi-colon, so this should be removed with second index '[0:-1]'.
    #         unifValue_file.close()
    #         os.system("rm unifValue")
    #         tempCoordFile = np.ones([nNodes, 1]) * float("inf")
    #         for j in np.arange(nNodes):
    #             tempCoordFile[j, 0] = unifValue
    #     if i == 0:
    #         coordList_tot = np.ones([nNodes, 4]) * float("inf")  # ID - X - Y - Z
    #         coordList_tot[:, 0] = np.arange(nNodes)  # CellID = numbers from 0 till (nNodes - 1)
    #     coordList_tot[:, (i + 1)] = tempCoordFile[:, 0]
    #
    # self.nNodes_tot = len(coordList_tot[:, 0])

    def Get_Point_IDs(self, boundary, dir):
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
            nFaces = int(line[:-2].split()[1])
            line = f.readline()
            startFace = int(line[:-2].split()[1])

        # Get number of points to keep a list of booleans
        prev_line = "("
        with open(f_p) as f:
            line = f.readline()
            while not "(" in line:
                prev_line = line
                line = f.readline()
            N_points = int(prev_line)
            points = np.zeros((N_points, 3))
            count = 0
            line = f.readline()
            while any(char.isdigit() for char in line):
                temp = line.split(" ")
                points[count, 0] = float(temp[0][1:])
                points[count, 1] = float(temp[1][:])
                points[count, 2] = float(temp[2][:-2])
                count += 1
                line = f.readline()

        # print(N_points)

        points_Bool = np.zeros((N_points, 1))
        boundary_Ind = []

        # Read in nodes file

        # Read in the list of faces and the nodes constituting those faces
        All_Fnodes = []
        with open(f_f) as f:
            line = f.readline()
            while not "(" in line:
                line = f.readline()
            line = f.readline()
            while any(char.isdigit() for char in line):
                list = line[2:-2].split()
                All_Fnodes.append(list)
                line = f.readline()

        # Extract the cell faces belonging to the boundary and create an ordered list of node id's
        face_centers = np.zeros((nFaces, 3))
        node_coords = []
        for nf in range(startFace, startFace + nFaces):
            face_center = np.zeros([1, 3])
            list = All_Fnodes[nf]
            for n in list:
                face_center += points[int(n), :]
                if points_Bool[int(n) - 1, 0] == 0:  # If not yet in list, add it to the list
                    boundary_Ind.append(int(n))
                    node_coords.append(points[int(n), :])
                    points_Bool[int(n) - 1, 0] = 1
            face_centers[nf - startFace, :] = face_center / float(len(list))
        boundary_Ind = np.array(boundary_Ind)
        node_coords = np.array(node_coords,dtype=float)

        return boundary_Ind, node_coords,face_centers
            
    
