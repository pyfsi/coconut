# Abaqus

This is the documentation for all Abaqus solver wrappers. Abaqus is a structural solver implementing the finite element method.
Currently this wrapper only supports FSI simulations, no other multi-physics problems. 
Subcycling within the structural solver is possible.

## Fluid-structure interaction with Abaqus
Abaqus (Dassault Systèmes) can be used to solve for the structural displacement/deformation in partitioned FSI-simulations. The FSI interface consist of a *surfaces* in the Abaqus model, where pressure and surface traction loads are applied, and corresponding node *sets*, where the resulting computed displacements are returned to the solver wrapper. The loads are applied in so-called load points (Gauss points, quadrature points), the displacements are exported in the elements' nodes. The input loads are collected in one or more `ModelParts` in the *input* `Interface`, the output nodes are collected in one or more `ModelParts` of the *output* `Interface`. Each `ModelPart` on the input `Interface` has a counterpart on the output `Interface`. More information about `ModelParts` and `Interface` can be found in the [data structure documentation](../../../data_structure/data_structure.md).

## Terminology
 - Main directory: Directory where the analysis is started.
 - Working directory: Subdirectory of the main directory in which Abaqus runs.
 - Source directory: Directory where the source files of the Abaqus solver wrapper are found: *`coupling_components/solver_wrappers/abaqus`*.
 - Extra directory: Subdirectory of the source directory with some files to assist with the setup.
 - Geometrical nodes: Nodes in Abaqus related to the geometry of the elements. At these nodes the displacement data is exported. 
 - Load points: Every element has load points. This is where the loads (input to Abaqus) are applied.
 - Time step: Time step from the viewpoint of the fluid-structure interaction, usually equal to the time-step of the flow solver and structural solver, although one of the solvers can deviate in the case of subcycling.
 - Increment: Time increment in the nomenclature of the Abaqus software. This is usually equal to the time step of the flow solver and overall coupled simulation, but in case of subcycling within the Abaqus solver, a time step can be subdivided in multiple increments.

## Environment
 - A working directory for Abaqus needs to be created within the main directory. Its **relative path to the main directory** should be specified in the JSON file. In the [CoCoNuT examples](../../../examples/examples.md) this folder is typically called *`CSM`*, but any name is allowed.
 - Abaqus needs a **host-file** called *`AbaqusHosts.txt`* in the working directory.
    - This host-file lists the machines on which Abaqus is allowed to run. One line per requested core, but excessive lines cause no harm.
    - The extra directory contains a script *`Example_makeHostFile.sh`* which can be used to generate a host file.
 - The **Abaqus license server** needs to be specified in the parametrized file `abaqus_v6.env` which is present in the `source directory`. For use at Ghent University no changes are required. 
 - The **Abaqus software should be available as well as compilers** to compile the user-subroutines (FORTRAN) and post-processing code (C++). Some compilers also require a license. 
     - The extra directory contains a file *`abaqus_setup`* which can be sourced in the terminal that will be used to run the simulation. This script loads the Abaqus module as well as the compilers (Ghent University system).

## Parameters
This section describes the parameter settings in the JSON file. A distinction is made between mandatory and optional parameters. It can be useful to have a look at a JSON file of one of the examples in the *`examples`* folder.

### Mandatory
parameter|type|description
---:|:---:|---
`arraysize`|int|Size specification for array in FORTRAN part of the code, to reserve sufficient memory. Should be large enough and depends on the number of load points in the structural model.
`cores`|int|Number of cores to be used by Abaqus.
`delta_t`|float|Size of the time step in Abaqus. Its value should be synchronized with the flow solver. This parameter is usually specified in a higher `Component` object in which case it is not mandatory.
`dimensions`|int|Dimensionality of the problem (2 or 3).
`surfaceIDs`|list|List with the names of the node sets associated with the geometrical nodes on the FSI-interface surfaces. <br><br> **Example:** ["NODESET_NAME_A", "NODESET_NAME_B"] <br> <br> **Important notes:**<br> &emsp;•	The sequence of these surfaces has to correspond with the integer specified in the corresponding load-surface (see the [input file section](#setting-up-a-case-abaqus-input-file-inp)).<br>&emsp;•	The names of these surfaces should be all UPPERCASE as Abaqus recasts these when opening the .inp file.
`interface_input`|list|Should contain a dictionary for each corresponding element in `surfaceIDs`, with a key `"model_part"` that provides the name of a `ModelPart` for Abaqus load points as value. The second key of the dictionary is `variables`. The list given as value specifies the input variables that should be included, chosen from *`data_structure/variables.py`*. Currently only `"pressure"` and `"traction"` are allowed (case-sensitive). An example can be found in [this part of the input file section](#input-related-settings-in-json-file).
`interface_output`|list|Similar to interface_input but for Abaqus geometrical nodes. In this case the `"variables"` key specifies the output variable, chosen from *`data_structure/variables.py`*. Currently only `"displacement"` is allowed (case-sensitive). An example can be found in [this part of the input file section](#output-related-settings-in-json-file).
`input_file`|str|Name of the Abaqus input file (.inp) provided by the user. <br> <br> **Example:** `"Base.inp"`
`mp_mode`|str|Determines how Abaqus is executed in parallel. Should be `"THREADS"` as `"MPI"` is currently not implemented.
`save_iterations`|int|Determines what files are kept by Abaqus. All files are saved, but files not corresponding to (i.e. of which the time step is not a multiple of) `save_iterations` are removed at the end of a time step. Important for restart options (also in correspondence with the save interval of the flow solver).
`timestep_start`|int|Time step to start from. Data should be available at this time step. For a new simulation this value will typically be 0. This parameter should be synchronized with the flow solver. This parameter is usually specified in a higher `Component` in which case it is not mandatory to specify. 
<nobr>`working_directory`</nobr>|str|Relative path to the directory in which Abaqus will be executed and where all structural information will be stored. Should be created before execution and contain a file *`AbaqusHosts.txt`*, see the [environment section](#environment).

`timestep_start` and `delta_t` are necessary parameters, but are usually defined in a higher `Component`. However, they can also be given directly as parameter of the solver wrapper (e.g. for standalone testing). If they are defined both in higher object and in the solver wrapper, then the former value is used and a warning is printed.

### Optional
parameter|type|description
---:|:---:|---
<nobr>`subcycling`</nobr>|boolean|`false`: [Default] Abaqus solves the requested time step using one increment. <br> `true`: Abaqus is allowed to solve the time step using multiple *increments*. This can be of use when Abaqus has convergence difficulties. For example cases where contact is involved often require small *increments*.
`initialInc`|float|Required when subcycling is enabled. Contains the size of the first time *increment* attempted by Abaqus.
`maxInc`|float|Required when subcycling is enabled. Contains the maximal time *increment* size allowed. This value should not be higher than `delta_t`.
`maxNumInc`|int|Required when subcycling is enabled. Contains the maximum number of *increments* that Abaqus is allowed to perform for one time step. 
`minInc`|float|Required when subcycling is enabled. Contains the minimal size allowed for a time *increment*.
`ramp`|boolean| Only used when subcycling is enabled in Abaqus. <br> `false`: Load is considered to be constant throughout the time step. <br>`true`: Load is applied in a ramped fashion throughout the time step. 

## Overview of operation
The solver wrapper consists of 5 files located in the source directory (with *`X`* denoting the Abaqus version, e.g. *`v614.py`*):

 - *`X.py`*: defines the `SolverWrapperAbaqusX`class. 
 - *`abaqus_v6.env`*: environment file setting the environment for the Abaqus solver.
 - *`GetOutput.cpp`*: Extracts the output (from Abaqus .odb files) and writes it to a file for each output`ModelPart`. Written in C++.
 - *`USR.f`*: An Abaqus user-subroutine that reads the loads from files (one for each input `ModelPart`) and applies them on the load points. Written in FORTRAN.
 - *`USRinit.f`*: An Abaqus user-subroutine that extract the coordinates of the load points and writes them to files (one for each input `ModelPart`) to initialize each input `ModelPart`. Written in FORTRAN.

### The `__init__` method

 During initialization of the `SolverWrapperAbaqusX` object, some parameters are substituted in *`abaqus_v6.env`*, *`GetOutput.cpp`*, *`USR.f`* and *`USRinit.f`* and these files are copied to the working directory. The C++ files and FORTRAN files are subsequently compiled. USRinit is ran to obtain the coordinates of the load points of each surfaces of which the name matches "MOVINGSURFACE**B**", where **B** should correspond to the index of the elements in `surfaceIDs` (0, 1, 2, ...). These coordinates are stored in `ModelParts` of which the name corresponds to the entries in `interface_input` (which also correspond to elements in `surfaceIDs`). GetOutput is ran to extract the coordinates of the geometrical nodes. These coordinates are added to `ModelParts` of which the names corresponds to entries of `interface_output` (also matched with `surfaceIDs`). The input `ModelParts` are added to an [`Interface`](../../../data_structure/data_structure.md) object taking care of the inputs (i.e. loads), the output `ModelParts` to another instance of `Interface`taking care of outputs (i.e. displacements).
 
#### Files written in the working directory during `__init__`
 
 In the file conventions *`A`* is the start time step (`timestep_start` in the JSON file) and *`B`* the index of the corresponding element of `surfaceIDs`.
 
 - The Abaqus input file (`input_file` in JSON file) is processed into a file *`CSM_TimeA.inp`* and `CSM_Restart.inp`, the latter taking care of all simulations (i.e. coupling iterations) but the first.
 - Upon running USRinit the load point coordinates of each surface of which the name matches "MOVINGSURFACEB" are written to *`CSM_TimeACpu0SurfaceBFaces.dat`* and *`CSM_TimeACpu0SurfaceBFacesBis.dat`*. When these are processed by the solver wrapper, also *`CSM_TimeASurfaceBElements.dat`* is created.
 - Upon running GetOutput the geometrical nodes are written to *`CSM_TimeASurfaceBNodes.dat`*. 
 
### The `solve_solution_step` method
 
 This method of `SolverWrapperAbaqusX` is called each coupling iteration with an `Interface` object (input_interface) containing loads, which are written to files that are read by the (compiled) *`USR.f`* during the invoked Abaqus simulation. The Abaqus software is started and shut down for each calculation, i.e. each coupling iteration. When the simulation has ran successfully (log-file *`AbaqusSolver.log`* is checked for errors), the outputs are read from Abaqus by GetOuput and written to a file. The file is read in Python and the output (displacements) are stored in the output_interface `Interface` object which is returned.
 
#### Files written in the working directory during `solve_solution_step`
 
In the file conventions *`A`* is the time step and *`B`* the index of the corresponding element of `surfaceIDs`.

 - Files written by Abaqus for allowing a restart (required every coupling iteration): *`CSM_TimeA.odb`*, *`CSM_TimeA.res`*, *`CSM_TimeA.mdl`*, *`CSM_TimeA.prt`*, *`CSM_TimeA.stt`*.
 - Output database file written by Abaqus called *`CSM_TimeA.odb`* and read by GetOutput (also needed for restart).
 - Output text file *`CSM_TimeASurfaceBOutput.dat`* containing displacements written by GetOutput and read by the solver wrapper.
 - Input text file *`CSM_TimeASurfaceBCpu0Input.dat`* containing the loads written by the solver wrapper and read by the USR.

## Setting up a case: Abaqus input file (.inp)
The Abaqus solver wrapper is configured to start from an input file which contains all necessary information for the calculation. This file should be located in the main directory. Its name should be specified in the JSON file via the parameter `input_file`. For the remainder of this section this file will be referred to as "base-file".

Creation of the base-file is not considered a part of the solver wrapper functionality as it is case-specific. However, in order for the solver wrapper to work, the base-file has to comply with certain general conditions. This section aims at informing the user about the requirements for the base-file.

### General
The base-file needs to be of the ".inp" type, this is an "input file for Abaqus". ".inp-files" are created via Abaqus by, after configuration, creating a "job" and requesting a "write input" for that job. These files can be opened in Abaqus by using "file > import > model".
The base-file has to contain all necessary information about the structural model, which includes:

 - Mesh defining the structure geometry and discretization.
    - Also the element type needs to be defined.
 - Material properties.
 - Boundary conditions.
 - Surfaces where external loads need to be applied.
    - Here *"Surface"* refers to nomenclature of the Abaqus software itself. The name can be found as such in the Abaqus working tree.
 - Node sets where displacement data will be extracted.
    - Here *"Set"* refers to nomenclature of the Abaqus software itself. The name can be found as such in the Abaqus working tree. Also element sets exist, but for CoCoNuT the demanded sets need to be a collection of geometrical nodes, hence "node set".
 - A *Step* definition, which contains solver settings. Currently the following type of analyses (it is advised to explicitly set them rather than leaving it to Abaqus to fill in a default) are supported:
    - Implicit dynamic, application quasi-static
    - Implicit dynamic, application moderate dissipation
    - Static general
 - Additional loads not dependent on the flow solver.
 
Abaqus models contain parts and those parts are used to create assemblies. The base-file should contain one assembly, which will then be used by the coupling. The assembly, thus, determines the position and orientation that will be used by the coupling software. Also the sets and surfaces required by CoCoNuT should be defined on the assembly level.
 
 Abaqus has a GUI as well as a Python 2 interface (which is also accessible) via the GUI. References to both the Python interface and GUI will be made below.
 
### Setup for Abaqus input (loads)
Per surface in the fluid-structure interface (where loads and displacements need to be exchanged) a "surface" should be created **in the assembly**. There are multiple possibilities to create these surfaces:

 - From the geometry: when the geometry has been defined in Abaqus itself, the geometry faces can easily be selected in the GUI. This method is often the most straightforward, but the Abaqus model should contain the geometry.
 - From the mesh: when the geometry is not available (this can for example be the case when a mesh has been imported), a surface can be defined by selecting multiple mesh faces. As a surface typically covers many mesh faces, it is useful there to select the regions "by angle", which uses the angle between mesh faces to determine whether adjacent faces should be selected. This way the surface selection can be extended until a sharp corner is met.
 - By converting a "node set" containing all the nodes on the surface and then calling the `SurfaceFromNodeSet` method which can be found in the `makeSurface.py` file in the extra directory.

The name of the surface has to be **MOVINGSURFACE followed by an integer**. The integer of the surface has to correspond with the index of that specific surface in the `surfaceIDs` parameter (index starts from 0). For example MOVINGSURFACE0 is associated with the first item in `surfaceIDs` and MOVINGSURFACE1 would be associated with the second item in `surfaceIDs`.

An example on the use of `SurfaceFromNodeSet` (via the Python console in Abaqus or a Python script for Abaqus):
```python
from makeSurface import SurfaceFromNodeSet
my_model = mdb.models['Model-1']
my_assembly = my_model.rootAssembly  
my_instance = my_assembly.instances['PART-1-1']
movingSurface0 = SurfaceFromNodeSet(my_model, my_instance, 'NODESET_NAME_A', 'MOVINGSURFACE0')
```

On these surfaces a "pressure load" and a "surface traction load" need to be specified with a "user-defined" distribution. Loads are assigned to a "step". A step is a part of the simulation to which an analysis type, algorithm settings and incrementation settings are assigned that do not change for the duration of the step. In a typical CoCoNuT case only a single step is defined. **Note that the name of the step used for the FSI has to be "Step-1" as this name is hardcoded in `GetOutput.cpp`**. After creation of the step the loads can be assigned. This can be done via the GUI or using Python commands available in Abaqus similar to the following:

```python
from step import *
step1 = my_model.ImplicitDynamicsStep(name='Step-1', previous='Initial', timePeriod=1, nlgeom=ON, maxNumInc=1, haftol=1, initialInc=1, minInc=1, maxInc=1, amplitude=RAMP, noStop=OFF, nohaf=ON, initialConditions=OFF, timeIncrementationMethod=FIXED, application=QUASI_STATIC)
step1.Restart(frequency = 99999, overlay = ON)
my_model.Pressure(name = 'DistributedPressure', createStepName = 'Step-1', distributionType = USER_DEFINED, field = '', magnitude = 1, region=movingSurface0)
my_model.SurfaceTraction(name = 'DistributedShear', createStepName = 'Step-1', region = movingSurface0, magnitude = 1, traction = GENERAL, directionVector = ((0,0,0), (1,0,0)), distributionType = USER_DEFINED)
```

The second command enables writing of restart files by Abaqus, which is required for running unsteady cases. This can be done from the GUI when the "Step" module is active in the viewport, by selecting "Output" in the top menu and subsequently "Restart Requests". *Frequency* should be put on *99999*, *overlay* *activated* (this spares disk space since only the last increment is kept) and _interval_ on _1_ (also see [this Abaqus documentation page](https://abaqus-docs.mit.edu/2017/English/SIMACAECAERefMap/simacae-t-simodbrestart.htm)).
Note that the step type "ImplicitDynamicsStep" is an example, it depends on the analysis procedure chosen to run the Abaqus part of the FSI simulation. For a steady simulation this could for example be "StaticStep":

```python
step1 = my_model.StaticStep(name='Step-1', previous='Initial', timePeriod=1.0, initialInc=1, minInc=1e-4, maxNumInc=10, nlgeom=ON, amplitude=RAMP)
```

Currently, the Abaqus wrapper automatically checks if the increments comply with the `delta_t`, `subcycling` and related settings and adjusts them accordingly or raises an error, but *only* for a *dynamic step* or a *static step with subcycling*.
**Attention:** Replacing the incrementation settings is done by looking up some keywords in the base-file (`*Step`, `*Dynamic`, `application`). This procedure fails when these keywords are not found. When using the GUI to create the base-file (.inp) and using the default settings for the step, often the keyword `application` is not written. It is hence advised not to use the default settings but use an application (quasi-static works for most cases and also moderate dissipation is allowed). If sub-cycling is not enabled, the maximal number of increments should be 1 (`inc=1`), otherwise an error is raised. This behavior may be changed in future versions of the code. The lines in the base-file (.inp) should look similar to this:

```
*Step, name=Step-1, nlgeom=YES, inc=1
*Dynamic,application=QUASI-STATIC,direct,nohaf,initial=NO
0.0001,0.0001,
```

The time step (0.0001) will in this case be replaced by settings found in the JSON file. More information can be found in [this Abaqus documentation page](https://abaqus-docs.mit.edu/2017/English/SIMACAEKEYRefMap/simakey-r-dynamic.htm).

#### Input-related settings in JSON file
Although it is imperative to use MOVINGSURFACE followed by an integer in the base-file to identify the surfaces *for input*, some lines in the JSON file (see example below) related to input contains names specified by the user. The names in the `interface_input["model_part"]` keys are matched by the elements (and their index) in `surfaceIDs`. 

```json
{
  "surfaceIDs": ["NODESET_NAME_A", "NODESET_NAME_B"],
  "interface_input": 
  [
    {
      "model_part": "NODESET_NAME_A_load_points",
      "variables": ["pressure", "traction"]
    },
    {
      "model_part": "NODESET_NAME_B_load_points",
      "variables": ["pressure", "traction"]
    }
  ]
}
```

### Setup for Abaqus output (displacements)
After creation of the step, Abaqus needs to be instructed about what to output at the end of a calculation. A "Field Output" has to be generated covering all locations involved in the fluid-structure interface. To do so one must create node sets in the assembly (if this had not been done before) with names corresponding to `surfaceIDs` and containing all structural nodes of the surfaces, then create a Field Output Request for at least the coordinates and the displacements. This Field Output Request can in the GUI be found as part of the model tree, but below also an example for the Python interface is given. A Field Output Request requests field output (as the name says) to be written to the output database file (.odb).

In the previous section an example was given of how a surface can be created from a node set, but the other way around is also possible, creating a node set from a surface (presuming that this surface was already created):

```python
my_assembly = my_model.rootAssembly
movingSurface0 = my_assembly.surfaces["MOVINGSURFACE0"]
outputSet = my_assembly.Set(name='NODESET_NAME_A', nodes=movingSurface0.nodes)
my_model.FieldOutputRequest(createStepName='Step-1', frequency=LAST_INCREMENT, name='F-Output-1', region=my_assembly.sets['NODESET_NAME_A'], variables=('COORD', 'U'))
```

Furthermore, it may be interesting (for post-processing and debugging) to preserve the default Abaqus output and therefore also configure a Field Output and History Output with PRESELECTED variables. This can be done via the GUI or using Python lines similar to the following:

```python
my_model.FieldOutputRequest(createStepName='Step-1', frequency=LAST_INCREMENT, name='F-Output-2', variables=PRESELECT)
my_model.HistoryOutputRequest(createStepName='Step-1', frequency=LAST_INCREMENT, name='H-Output-1', variables=PRESELECT)
```

#### Output-related settings in JSON file
The values of the `interface_output["model_part"]` keys are matched by the elements (and their index) in `surfaceIDs`. These values are internally used in CoCoNuT to distinguish the different `ModelParts`. The elements of `surfaceIDs` should match the names of the node sets defined in Abaqus.

```json
{
  "surfaceIDs": ["NODESET_NAME_A", "NODESET_NAME_B"],
  "interface_output": 
  [
    {
      "model_part": "NODESET_NAME_A_nodes",
      "variables": ["displacement"]
    },
    {
      "model_part": "NODESET_NAME_B_nodes",
      "variables": ["displacement"]
    }
  ]
}
```

### Note about choosing `ModelParts`
The created "surfaces" and "node sets" for load input and displacement output respectively, correspond to `ModelParts` in the CoCoNuT code, a representation of the data used for the coupling. It is strongly advised to sub-divide to fluid-structure interaction interface intelligently, depending on the geometry. As a rule of thumb it can be said that a surfaces at two sides of a sharp corner should be assigned to a different `ModelPart`. As the interpolation is based on shortest distance, issues can arise at sharp corners. Those are avoided by having different `ModelParts` at each side of the corner. Another reason to do this is because the code cannot handle elements with two or more faces being part of the same `ModelPart`. This situation would occur if the surface contains corners. An example is an airfoil where the suction side and pressure side belong to the same `ModelPart`: elements at the trailing edge will have (a) face(s) at both the pressure side and suction side. Even when the code would allow this, interpolation mistakes become likely, as a geometrical node or load point on the suction side could have a nearest neighbour on the pressure side, causing that the wrong data is used for interpolation.

## Log files

A general event log of the procedure can be found in the working directory, in a file named *`AbaqusSolver.log`*. For more detailed information on a certain time step, the .msg file written by Abaqus can be consulted. In CoCoNuT these are structured as follows: *`CSM_TimeA.msg`*, *`A`* being the time step. Typically multiple coupling iterations are done within each time step, so these .msg-files get overwritten by each new coupling iteration in the same time step.