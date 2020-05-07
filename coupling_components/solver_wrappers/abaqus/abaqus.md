# SolverWrapperAbaqus

This is the documentation for all Abaqus SolverWrappers.
Currently only FSI simulations are supported, no other multiphysics problems. 
Subcycling within the structural solver is possible.

## Terminology
`Main directory`: Directory where the Python code will be executed  
`Source directory`(for wrapper): <span style="color:black">.../coconut/coupling_components/solver_wrappers/abaqus/</span>  
`Extra directory`: Subdirectory of the source directory with some files to assist with the setup  
`Working directory`: Subdirectory of the main directory where Abaqus will be running   
`Geometrical nodes`: Nodes in Abaqus related to the geometry of the elements. At these nodes the displacement data is exported.   
`Loadpoints`: Every element has load points. This is where the loads (input to Abaqus) are applied.  

## Environment
 - A working directory for Abaqus needs to be created within the main directory. Its **relative path to the main directory** should be specified in the .json file.
 - Abaqus needs an `AbaqusHosts.txt` file in the working directory.
 - This host-file lists the machines on which Abaqus is allowed to run. One line per requested core, but excessive lines cause no harm.
 - The `extra directory` contains an `Example_makeHostFile.sh` which was used the generate a host file for the local system.
 - The Abaqus license server needs to be specified in the parametrized file `abaqus_v6.env` which is present in the `source directory` 
 - The `extra directory` contains a file `abaqus_setup` which should be sourced in the terminal that will be used to run the simulation. This script loads the Abaqus module as well as the compilers.


## Parameters 
This section describes the Parameters in the JSON file, listed in alphabetical order. A distinction is made between mandatory and optional parameters.
### Mandatory
parameter|type|description
---:|:---:|---
<nobr>`arraysize`</nobr>|integer|Size specification for array in FORTRAN part of the code. Should be large enough and depends on the number of load points in the structural model.
<nobr>`cores`</nobr>|integer|Number of cores to be used by Abaqus.
<nobr>`delta_t`</nobr>|float|Size of the time step in Abaqus (Should be synchronized with the flow solver). This parameter is usually specified in a higher `Component`.
<nobr>`dimensions`</nobr>|integer|Dimensionality of the problem (2 or 3).
<nobr>`interface_input`</nobr>|dict|Should contain `surfaces` keys. Keys are names of ModelParts for Abaqus load points. Each name (key) must be the contain an entry from `surfaceIDs`. The list that follows specifies the historical variables that should be included in this modelpart. (Note the comma) <br> <br> <b>Example:</b> <br> &emsp;"NODESETA_load_points": ["PRESSURE", "TRACTION"],<br>&emsp;"NODESETB_load_points": ["PRESSURE", "TRACTION"]
<nobr>`interface_output`</nobr>|dict|Similar to interface_input but for Abaqus load points.
<nobr>`input_file`</nobr>|string|Name of the input file provided by the user. <br> <b>Example:</b> “Base.inp”
<nobr>`mp_mode`</nobr>|string|Determines how Abaqus is executed in parallel. Should be “THREADS” as “MPI”  is currently not implemented
<nobr>`save_iterations`</nobr>|integer|Determines what files are kept by Abaqus. All files are saved, but files not corresponding to the save_iterations are removed at the end of a time step. Important for restart options (also in correspondence with the save interval of the flow solver).
<nobr>`surfaceIDs`</nobr>|list of strings|Comma-separated list with the names of the node sets associated with the geometrical nodes on the interface surfaces. <br><br> <b>Example:</b>  [“NODESETA”, “NODESETB”] <br> <br> <b>Important notes:</b><br> &emsp;•	The sequence of these surfaces has to correspond with the integer specified in the corresponding load-surface (See the Input file section)<br>&emsp;•	The names of these surfaces should be all UPPERCASE as Abaqus recasts these when opening the .inp file.
<nobr>`surfaces`</nobr>|integer|Number of surfaces on the FSI interface. Should correspond with the number of elements in `surfaceIDs`.
<nobr>`timestep_start`</nobr>|integer|Time step to start from [data should be available at this time step. For a new simulation this value will typically be 0] (Should be synchronized with the flow solver). This parameter is usually specified in a higher `Component`. 
<nobr>`working_directory`</nobr>|string|Relative path to the directory in which Abaqus will be executed and where all structural information will be stored. <br> Should be created before execution and contain an `AbaqusHosts.txt` file.

`timestep_start` and `delta_t` are necessary parameters, but are usually defined in a higher `Component`. However, they can also be given directly as parameter of the solver wrapper (e.g. for standalone testing). If they are defined both in higher object and in the solver wrapper, then the former value is used and a warning is printed.

If different parameters are used with different Abaqus versions, this should be specified both in this section and in the version specific documentation section.

### Optional
parameter|type|description
---:|:---:|---
<nobr>`initialInc`</nobr>|float|Required when subcycling is enabled. Contains the size of the first substep attempted Abaqus
<nobr>`maxInc`</nobr>|float|Required when subcycling is enabled. Contains the maximal size allowed for a substep
<nobr>`maxNumInc`</nobr>|integer|Required when subcycling is enabled. Contains the maximum number of increments that Abaqus is allowed to perform for one time step 
<nobr>`minInc`</nobr>|float|Required when subcycling is enabled. Contains the minimal size allowed for a substep
<nobr>`ramp`</nobr>|integer| Only important when subcycling is enabled in Abaqus. <br> <b>0</b>: Load is considered to be constant throughout the time step <br><b>1</b>: [Default] Load is applied in a ramped fashion throughout the time step 
<nobr>`subcycling`</nobr>|integer|<b>0</b>: [Default] Abaqus solves the requested time step using one increment <br> <b>1</b>: Abaqus is allowed to solve the time step using multiple increments, often required when contact is involved 

## Input file
The Abaqus solver wrapper is configured to start from an input file which contains all necessary information for the calculation. This file should be located in the `main directory`. Its name should be specified in the .json file via the parameter `input_file`. For the remainder of this section this file will be referred to as “base-file”
<br><br>Creation of the base-file is not considered a part of the solver wrapper functionality as it is case-specific. However, in order for the solver wrapper to work, the base-file has to comply with certain general conditions. This section aims at informing the user about the requirements for the base-file.

### General
The base-file needs to be of the “.inp” type, this is an “input file for Abaqus”. “.inp-files” are created via Abaqus by, after configuration, creating a “job” and requesting a “write input” for that job. These files can be opened in Abaqus by using “file > import > model”.
The base-file has to contain all necessary information about the structural model, which includes:  

 - Geometry
 - Mesh
 - Material properties
 - Boundary conditions
 - Surfaces where external loads need to be applied
 - Additional loads not dependent on the flow solver
 
Abaqus models contain parts and those parts are used to create assemblies. The base-file should contain one assembly, which will then be used by the coupling. The assembly, thus, determines the position and orientation that will be used by the coupling software.
 
### Setup for Abaqus input (loads)
Per surface in the fluid-structure interface (where loads and displacements need to be exchanged) a “surface” should be created in the assembly. 
These surfaces can for example be created from geometry or by creating a “node set” containing all the nodes on the surface and then calling the “SurfaceFromNodeSet” function which can be found in the makeSurface.py file in the `Extra directory`. The name of the surface has to be MOVINGSURFACE followed by an integer. The integer of the surface has to correspond with the index of that specific surface in the `surfaceIDs` parameter (index starts from 0). 
For example MOVINGSURFACE0 is associated with the first item in `surfaceIDs` and MOVINGSURFACE1 would be associated with the second item in `surfaceIDs`. 
An example on the use of SurfaceFromNodeSet (via the Python console in Abaqus or a python script for Abaqus):  
 
```python
from makeSurface import *
my_model=mdb.models['Model-1']
my_assembly=my_model.rootAssembly  
my_instance=my_assembly.instances['PART-1-1']
movingSurface0 = SurfaceFromNodeSet(my_model, my_instance, 'NAME_OF_THE_NODESET', 'MOVINGSURFACE0')
```  

On these surfaces a pressure load and a traction load need to be specified with a user-defined distribution. Loads are assigned to a “step”. After creation of the step the loads can be assigned. This can be done via the GUI or using python commands available in Abaqus similar to the following:  

```python
from step import *
step1 = my_model.ImplicitDynamicsStep(name='Step-1', previous='Initial', timePeriod=1, nlgeom=ON, maxNumInc=1, haftol=1, initialInc=1, minInc=1, maxInc=1, amplitude=RAMP, noStop=OFF, nohaf=ON, initialConditions=OFF, timeIncrementationMethod=FIXED, application=QUASI_STATIC)
my_model.Pressure(name = 'DistributedPressure', createStepName = 'Step-1', distributionType = USER_DEFINED, field = '', magnitude = 1, region=movingSurface0)
my_model.SurfaceTraction(name = 'DistributedShear', createStepName = 'Step-1', region = movingSurface0, magnitude = 1, traction = GENERAL, directionVector = ((0,0,0), (1,0,0)), distributionType = USER_DEFINED)
```

### Setup for Abaqus output (displacements)
After creation of the step Abaqus needs to be instructed about what to output at the end of a calculation. A fieldOutput 
has to be generated for each surface involved in the fluid-structure interface. 
To do so it is best to create node sets in the assembly containing all structural nodes of the surfaces (if this hadn’t been done before) 
and to create a fieldOutput per surface containing at least the coordinates and the displacements. 
Furthermore, it is interesting to output maintain the default Abaqus output and therefore also configure a fieldOutput 
and historyOutput with PRESELECTED variables. 
This can be done via the GUI or using python lines similar to the following:  

```python   
my_model.FieldOutputRequest(createStepName='Step-1', frequency=LAST_INCREMENT, name='F-Output-1', region=tubeAssembly.sets['NAME_OF_THE_NODE_SET'], variables=('COORD', 'U'))
my_model.FieldOutputRequest(createStepName='Step-1', frequency=LAST_INCREMENT, name='F-Output-2', variables=PRESELECT)
my_model.HistoryOutputRequest(createStepName='Step-1', frequency=LAST_INCREMENT, name='H-Output-1', variables=PRESELECT)
```  
