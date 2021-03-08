# Abaqus

This is the documentation for all Abaqus solver-wrappers. Abaqus is a structural solver implementing the finite element method.
Currently this wrapper only supports FSI simulations, no other multiphysics problems. 
Subcycling within the structural solver is possible.

## Terminology
-   Main directory: Directory where the analysis is started.
-   Working directory: Subdirectory of the main directory in which Abaqus runs.
-   Source directory: Directory where the source files of the Abaqus solver-wrapper are found: *`$COCO/coconut/coupling_components/solver_wrappers/abaqus`*.
-   Extra directory: Subdirectory of the source directory with some files to assist with the setup.
-   Geometrical nodes: Nodes in Abaqus related to the geometry of the elements. At these nodes the displacement data is exported.   
-   Load points: Every element has load points. This is where the loads (input to Abaqus) are applied.
-   Time step: Time step from the viewpoint of the fluid-structure interaction, usually equal to the time-step of the flow solver and structural solver, although one of the solvers can deviate in the case of subcycling.
-   Increment: Time increment in the nomenclature of the Abaqus software. This is usually equal to the time step of the flow solver and overall coupled simulation, but in case of subcycling within the Abaqus solver, a time step can be subdivided in multiple increments.

## Environment
 - A working directory for Abaqus needs to be created within the main directory. Its **relative path to the main directory** should be specified in the JSON file.  In the [CoCoNuT examples](../../../examples/examples.md) this folder is typically called *`CSM`*, but any name is allowed.
 - Abaqus needs a **host-file** called *`AbaqusHosts.txt`* in the working directory.
     - This host-file lists the machines on which Abaqus is allowed to run. One line per requested core, but excessive lines cause no harm.
     - The extra directory contains a script *`Example_makeHostFile.sh`* which can be used to generate a host file.
 - The **Abaqus license server** needs to be specified in the parametrized file `abaqus_v6.env` which is present in the `source directory`. For use at Ghent University no changes are required. 
 - The **Abaqus software should be available as well as compilers** to compile the user-subroutines (FORTRAN) and post-processing code (C++). Some compilers also require a license. 
     - The extra directory contains a file *`abaqus_setup`* which can be sourced in the terminal that will be used to run the simulation. This script loads the Abaqus module as well as the compilers (Ghent University system).

## Parameters
This section describes the parameter settings in the JSON file. A distinction is made between mandatory and optional parameters. It can be useful to have a look at a JSON file of one of the examples in *`$COCO/coconut/examples`*.

### Mandatory
parameter|type|description
---:|:---:|---
`arraysize`|int|Size specification for array in FORTRAN part of the code, to reserve sufficient memory. Should be large enough and depends on the number of load points in the structural model.
`cores`|int|Number of cores to be used by Abaqus.
`delta_t`|float|Size of the time step in Abaqus. Its value should be synchronized with the flow solver. This parameter is usually specified in a higher `Component` object in which case it is not mandatory.
`dimensions`|int|Dimensionality of the problem (2 or 3).
`surfaceIDs`|list|List with the names of the node sets associated with the geometrical nodes on the interface surfaces. <br><br> <b>Example:</b>  [“NODESETA”, “NODESETB”] <br> <br> <b>Important notes:</b><br> &emsp;•	The sequence of these surfaces has to correspond with the integer specified in the corresponding load-surface (see the Input file section).<br>&emsp;•	The names of these surfaces should be all UPPERCASE as Abaqus recasts these when opening the .inp file.
`surfaces`|int|Number of surfaces on the FSI interface. Should correspond with the number of elements in `surfaceIDs`.
`interface_input`|list|Should contain `surfaces` elements. Each element is a dictionary with a key `"model_part"`, containing the name of a `ModelPart` for Abaqus load points. Each name must contain an entry from `surfaceIDs`. The second key of the dictionary is `variables`. The list given as value specifies the input variables that should be included. Currently only `"pressure"` and `"traction"` are allowed (case-sensitive). 
`interface_output`|list|Similar to interface_input but for Abaqus geometrical nodes. In this case the `"variables"` key specifies the output variable. Currently only `"displacement"` is allowed (case-sensitive).
`input_file`|str|Name of the Abaqus input file (.inp) provided by the user. <br> <b>Example:</b> `"Base.inp"`
`mp_mode`|str|Determines how Abaqus is executed in parallel. Should be `"THREADS"` as `"MPI"`  is currently not implemented.
`save_iterations`|int|Determines what files are kept by Abaqus. All files are saved, but files not corresponding to (i.e. of which the time step is not a multiple of) `save_iterations` are removed at the end of a time step. Important for restart options (also in correspondence with the save interval of the flow solver).
`timestep_start`|int|Time step to start from. Data should be available at this time step. For a new simulation this value will typically be 0. This parameter should be synchronized with the flow solver. This parameter is usually specified in a higher `Component` in which case it is not mandatory to specify. 
`working_directory`|str|Relative path to the directory in which Abaqus will be executed and where all structural information will be stored. <br> Should be created before execution and contain a file *`AbaqusHosts.txt`*.

`timestep_start` and `delta_t` are necessary parameters, but are usually defined in a higher `Component`. However, they can also be given directly as parameter of the solver wrapper (e.g. for standalone testing). If they are defined both in higher object and in the solver wrapper, then the former value is used and a warning is printed.

### Optional
parameter|type|description
---:|:---:|---
`subcycling`|boolean|`false`: [Default] Abaqus solves the requested time step using one increment. <br> `true`: Abaqus is allowed to solve the time step using multiple increments, often required when contact is involved. 
`initialInc`|float|Required when subcycling is enabled. Contains the size of the first time *increment* attempted by Abaqus.
`maxInc`|float|Required when subcycling is enabled. Contains the maximal time *increment* size allowed. This value should not be higher than `delta_t`.
`maxNumInc`|int|Required when subcycling is enabled. Contains the maximum number of *increments* that Abaqus is allowed to perform for one time step. 
`minInc`|float|Required when subcycling is enabled. Contains the minimal size allowed for a time *increment*..
`ramp`|boolean| Only used when subcycling is enabled in Abaqus. <br> `false`: Load is considered to be constant throughout the time step. <br>`true`: Load is applied in a ramped fashion throughout the time step. 

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
 - Node sets where displacement data will be extracted
 - Additional loads not dependent on the flow solver
 
Abaqus models contain parts and those parts are used to create assemblies. The base-file should contain one assembly, which will then be used by the coupling. The assembly, thus, determines the position and orientation that will be used by the coupling software.
 
### Setup for Abaqus input (loads)
Per surface in the fluid-structure interface (where loads and displacements need to be exchanged) a “surface” should be created **in the assembly**. 
These surfaces can for example be created from geometry or by creating a “node set” containing all the nodes on the surface and then calling the `SurfaceFromNodeSet` method which can be found in the `makeSurface.py` file in the `Extra directory`. The name of the surface has to be **MOVINGSURFACE followed by an integer**. The integer of the surface has to correspond with the index of that specific surface in the `surfaceIDs` parameter (index starts from 0). 
For example MOVINGSURFACE0 is associated with the first item in `surfaceIDs` and MOVINGSURFACE1 would be associated with the second item in `surfaceIDs`. 
An example on the use of SurfaceFromNodeSet (via the Python console in Abaqus or a python script for Abaqus):  
 
```python
from makeSurface import *
my_model=mdb.models['Model-1']
my_assembly=my_model.rootAssembly  
my_instance=my_assembly.instances['PART-1-1']
movingSurface0 = SurfaceFromNodeSet(my_model, my_instance, 'NAME_OF_THE_NODESET', 'MOVINGSURFACE0')
```  

On these surfaces a pressure load and a traction load need to be specified with a user-defined distribution. Loads are assigned to a “step”. A step is a part of the simulation to which an analysis type, algorithm settings and incrementation settings are assigned that do not change for the duration of the step. In a typical CoCoNuT case only a single step is defined. After creation of the step the loads can be assigned. This can be done via the GUI or using python commands available in Abaqus similar to the following:  

```python
from step import *
step1 = my_model.ImplicitDynamicsStep(name='Step-1', previous='Initial', timePeriod=1, nlgeom=ON, maxNumInc=1, haftol=1, initialInc=1, minInc=1, maxInc=1, amplitude=RAMP, noStop=OFF, nohaf=ON, initialConditions=OFF, timeIncrementationMethod=FIXED, application=QUASI_STATIC)
step1.Restart(frequency = 99999, overlay = ON)
my_model.Pressure(name = 'DistributedPressure', createStepName = 'Step-1', distributionType = USER_DEFINED, field = '', magnitude = 1, region=movingSurface0)
my_model.SurfaceTraction(name = 'DistributedShear', createStepName = 'Step-1', region = movingSurface0, magnitude = 1, traction = GENERAL, directionVector = ((0,0,0), (1,0,0)), distributionType = USER_DEFINED)
```

The second command enables writing of restart files by Abaqus, which is required for running unsteady cases. This can be done from the GUI when the "step" module is active, by selecting "Output" in the top menu and subsequently "Restart Requests". _Frequency_ should be put on _99999_, _overlay_ _activated_ (this saves data since only the last increment is kept) and _interval_ on _1_ (also see [this Abaqus documentation page](https://abaqus-docs.mit.edu/2017/English/SIMACAECAERefMap/simacae-t-simodbrestart.htm)).  
Note that the step type "ImplicitDynamicStep" is an example, it depends on the analysis procedure chosen to run the Abaqus part of the FSI simulation. For a steady simulation this could for example be "StaticStep":

```python
step1 = my_model.StaticStep(name='Step-1', previous='Initial', timePeriod=1.0, initialInc=1, minInc=1e-4, maxNumInc=10, nlgeom=ON, amplitude=RAMP)
```

Currently, the Abaqus wrapper automatically checks if the increments comply with the `delta_t` and `subcycling` settings and adjusts them accordingly or raises an error, but _only_ for a _dynamic step_ or a _static step with subcycling_.  
**Attention:** Replacing the incrementation settings is done by looking up some keywords in the base file (`*Step`, `*Dynamic`, `application`). This procedure fails when these keywords are not found. When using the GUI to create the base file (.inp) and using the default settings for the step, often the keyword `application` is not written. It is hence advised not to use the default settings but use an application (quasi-static works for most cases and also moderate dissipation is allowed). If sub-cycling is not enabled, the maximal number of increments should be 1 (`inc=1`), otherwise an error is raised. This behavior may be changed in future versions of the code. The lines in the input file should look similar to this:

```
*Step, name=Step-1, nlgeom=YES, inc=1
*Dynamic,application=QUASI-STATIC,direct,nohaf,initial=NO
0.0001,0.0001,
```

The time step (0.0001) will in this case be replaced by settings found in the json-file. More information can be found in [this Abaqus documentation page](https://abaqus-docs.mit.edu/2017/English/SIMACAEKEYRefMap/simakey-r-dynamic.htm).

### Setup for Abaqus output (displacements)
After creation of the step Abaqus needs to be instructed about what to output at the end of a calculation. A fieldOutput has to be generated covering all locations involved in the fluid-structure interface. 
To do so it is best to create node sets in the assembly containing all structural nodes of the surfaces (if this had not been done before) and to create a fieldOutput per surface containing at least the coordinates and the displacements. 

In the previous section an example was given of how a surface can be created from a node set, but the other way around is also possible, creating a node set from a surface (presuming that this surface was already created):

```python
my_model = mdb.models['Model-1']
my_assembly = my_model.rootAssembly  
movingSurface0 = my_assembly.surfaces["MOVINGSURFACE0"]
outputSet = my_assembly.Set(name='NAME_OF_THE_NODESET', nodes=movingSurface0.nodes)
my_model.FieldOutputRequest(createStepName='Step-1', frequency=LAST_INCREMENT, name='F-Output-1', region=my_assembly.sets['NAME_OF_THE_NODE_SET'], variables=('COORD', 'U'))
```

Furthermore, it is interesting (for post-processing and debugging) to preserve the default Abaqus output and therefore
also configure a fieldOutput and historyOutput with PRESELECTED variables. 
This can be done via the GUI or using python lines similar to the following:  

```python   
my_model.FieldOutputRequest(createStepName='Step-1', frequency=LAST_INCREMENT, name='F-Output-2', variables=PRESELECT)
my_model.HistoryOutputRequest(createStepName='Step-1', frequency=LAST_INCREMENT, name='H-Output-1', variables=PRESELECT)
```  

## ModelParts
The created "surfaces" and "node sets" for load input and displacement output respectively, correspond to ModelParts in the CoCoNuT code, a representation of the data used for the coupling. It is strongly advised to sub-divide to fluid-structure interaction interface intelligently, depending on the geometry. As a rule of thumb it can be said that a surfaces at two sides of a sharp corner should be assigned to a different ModelPart. As the interpolation is based on shortest distance, issues can arise at sharp corners. Those are avoided by having different ModelParts at each side of the corner.  
Another reason to do this is because the code cannot handle elements with two or more faces being part of the same ModelPart. This situation would occur if the surface contains corners. An example is an airfoil where the suction side and pressure side belong to the same ModelPart: elements at the trailing edge will have (a) face(s) at both the pressure side and suction side. Even when the code would allow this, interpolation mistakes become likely, as a `Node` on the suction side could have a nearest neighbour on the pressure side, causing that the wrong data is used for interpolation.
