# AbaqusCSE

This is the documentation for the new version of the Abaqus solver wrapper, which uses the Abaqus Co-Simulation Engine (CSE) avoiding the need to start and stop Abaqus every coupling iteration.
Refer to [Abaqus](../abaqus.md) for the first version of the Abaqus solver wrapper.

Abaqus is a structural solver implementing the finite element method.
Currently, this wrapper only supports FSI simulations, no other multi-physics problems.
Subcycling has not been fully tested but is expected to work within the structural solver.

!!! warning "For two-dimensional cases, the traction forces are not taken into account"

## Fluid-structure interaction with Abaqus
Abaqus (Dassault Systèmes) can be used to solve for the structural displacement/deformation in partitioned FSI-simulations.
The FSI interface consist of a *surfaces* in the Abaqus model, where pressure and traction loads are applied.
The loads are applied in the element centers, the displacements are exported in the element nodes.
The input loads are collected in one or more `ModelParts` in the *input* `Interface`,
the output nodes are collected in one or more `ModelParts` of the *output* `Interface`.
Each `ModelPart` on the input `Interface` has a counterpart on the output `Interface`.
More information about `ModelParts` and `Interface` can be found in the [data structure documentation](../../../data_structure/data_structure.md).

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
 - The **Abaqus software should be available as well as compilers** to compile the user-subroutines (FORTRAN) and post-processing code (C++). Some compilers also require a license. 
 - If the **Abaqus license server** needs to be specified explicitly, it is advised to do this in the [solver modules file](../../../README.md#checking-the-solver-modules).

## Parameters
This section describes the parameter settings in the JSON file.
It can be useful to have a look at a JSON file of one of the examples in the *`examples`* folder.

|                                         parameter | type | description                                                                                                                                                                                                                                                                                                                                                                                                                |
|--------------------------------------------------:|:----:|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|                                           `cores` | int  | Number of cores to be used by Abaqus.                                                                                                                                                                                                                                                                                                                                                                                      |
|                                           `debug` | bool | (optional) Default: `false`. Additional files are created or retained that are helpful for debugging, such as Abaqus send-and-receive file (Abaqus.SRE) or the input and output data of the solver for every coupling iterations.                                                                                                                                                                                          |
|                                      `dimensions` | int  | Dimensionality of the problem (2 or 3).                                                                                                                                                                                                                                                                                                                                                                                    |
| <nobr>`disable_modification_of_input_file`</nobr> | bool | (optional) Default: `false`. If `true`, the input file is not modified by CoCoNuT but simply copied to a file named Abaqus.inp. The user is responsible to correctly define the co-simulation settings, the time step, the duration, the output and the restart behaviour.                                                                                                                                                 |
|                                 `interface_input` | list | Should contain a dictionary for each input `ModelPart` to be created, having a key `model_part` and a key `variables`. The order of the list should correspond to the order of the Abaqus *surfaces* provided in the parameter `surfaces`, as well as to the order of the model parts in the other solver. The value of the key `variables` should be a list of the variables `pressure` and `traction`.                   |
|                                `interface_output` | list | Similar to `interface_input` but contains the output `ModelParts` for Abaqus geometrical nodes. In this case the `variables` key specifies the output variable, chosen from *`data_structure/variables.py`*. Here, the value of the key `variables` should be a `[displacement]`. Again, the order should correspond to the `surfaces` as well as to the `Interface` definitions of the solver to which Abaqus is coupled. |
|                                      `input_file` | str  | Name of the Abaqus input file (.inp) provided by the user. <br> Example: *`case.inp`*                                                                                                                                                                                                                                                                                                                                      |
|                                            `port` | int  | (optional) If not provided, the operating system will choose a free port. This port is used for the communication between the AbaqusWrapper executable and the CSE.                                                                                                                                                                                                                                                        |
|                                    `save_results` | int  | (optional) Default: `1`. This will modify the frequency of storing fields or history in the *`.odb`* file. If no output definition is present in the input file, one is created for the *PRESELECT* variables.                                                                                                                                                                                                             |
|                                        `surfaces` | list | List of the surface names defined in Abaqus. Note that the order has to be the same as those of the interface definitions. Take the Abaqus naming conventions into account: surfaces names are capitiliazed unless specified within double quotes, and the name is often prepended for example with ASSEMBLY_ (“Assembly.Part” is translated to “Assembly_Part”).                                                          |
|                               `working_directory` | str  | Relative path to the directory in which Abaqus will be executed and where all structural information will be stored. Should be created before execution and contain the input file.                                                                                                                                                                                                                                        |

The following parameters are usually defined in the top-level settings of the JSON file, but they can also be given directly as parameter of the solver wrapper (e.g. for standalone testing).
If they are defined in both locations, a warning is printed and the top-level value is used.

|                          parameter | type  | description                                                                                                                                              |
|-----------------------------------:|:-----:|----------------------------------------------------------------------------------------------------------------------------------------------------------|
|                          `delta_t` | float | Size of the time step in Abaqus.                                                                                                                         |
| <nobr>`number_of_timesteps`</nobr> |  int  | The amount of time steps to run the calculation. For a steady calculation, the value should be `1`.                                                      |
|                     `save_restart` |  int  | Indicates the time step interval at which files for restart have to be saved. A minus sign indicates only the files from the last interval are retained. |
|                   `timestep_start` |  int  | Time step to start from. Data should be available at this time step. For a new simulation this value will typically be `0`.                              |


## Overview of operation
The solver wrapper consists of 6 types of files located in the source directory (with *`X`* denoting the Abaqus version, e.g. *`v2022.py`*):

 - *`abaqus.py`*: Contains the base class `SolverWrapperAbaqus`.
 - *`X.py`*: Defines the `SolverWrapperAbaqusX`class, which inherits from the base class. Some version specific parameters might be overwritten in these subclasses. 
 - *`abaqus_v6.env`*: Environment file setting the environment for the Abaqus solver.


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
 - Surfaces where external loads need to be applied (one surface per `ModelPart`).
    - Here *"Surface"* refers to nomenclature of the Abaqus software itself. The name can be found as such in the Abaqus working tree.
 - Per surface a pressure load and traction load should be defined (see [below](#setup-for-abaqus-input-loads)).
 - Node sets where displacement data will be extracted.
    - Here *"Set"* refers to nomenclature of the Abaqus software itself. The name can be found as such in the Abaqus working tree. Also element sets exist, but for CoCoNuT the demanded sets need to be a collection of geometrical nodes, hence "node set".
- A Field Output Request requesting output at the node sets (see [below](#setup-for-abaqus-output-displacements)).
 - A *Step* definition, which contains solver settings. Currently, the following type of analyses are supported (it is advised to explicitly set them based on this documentation rather than leaving it to Abaqus to fill in a default):
    - Implicit dynamic, application quasi-static
    - Implicit dynamic, application moderate dissipation
    - Implicit dynamic, application transient fidelity
    - Static general
 - Additional loads not dependent on the flow solver.
 
Abaqus models contain parts and those parts are used to create assemblies. The base-file should contain one assembly, which will then be used by the coupling. The assembly, thus, determines the position and orientation that will be used by the coupling software. Also the sets and surfaces required by CoCoNuT should be defined on the assembly level.
 
 Abaqus has a GUI as well as a Python 2 interface (which is also accessible) via the GUI. References to both the Python interface and GUI will be made below.
 
### Setup for Abaqus input (loads)
Per surface in the fluid-structure interface (where loads and displacements need to be exchanged) a "surface" should be created **in the assembly**. There are multiple possibilities to create these surfaces:

 - From the geometry: when the geometry has been defined in Abaqus itself, the geometry faces can easily be selected in the GUI. This method is often the most straightforward, but the Abaqus model should contain the geometry.
 - From the mesh: when the geometry is not available (this can for example be the case when a mesh has been imported), a surface can be defined by selecting multiple mesh faces. As a surface typically covers many mesh faces, it is useful there to select the regions "by angle", which uses the angle between mesh faces to determine whether adjacent faces should be selected. This way the surface selection can be extended until a sharp corner is met.
 - By converting a "node set" containing all the nodes on the surface and then calling the `SurfaceFromNodeSet` method which can be found in the `make_surface.py` file in the extra directory. It is advised to copy this file to the working directory.

An example on the use of `SurfaceFromNodeSet` (via the Python console in Abaqus or a Python script for Abaqus):
```python
from make_surface import SurfaceFromNodeSet
my_model = mdb.models['Model-1']
my_assembly = my_model.rootAssembly  
my_instance = my_assembly.instances['PART-1-1']
inputSurfaceA = SurfaceFromNodeSet(my_assembly, my_instance, 'NODESET_NAME_A', 'SURFACE_NAME_A')
```

A step is a part of the simulation to which an analysis type, algorithm settings and incrementation settings are assigned that do not change for the duration of the step.
In a typical CoCoNuT case only a single step is defined.
After creation of the step the loads can be assigned.
It is required to have the `initialConditons=OFF` part when using `application=TRANSIENT FIDELITY` (beware that when application is not specified, transient fidelity is the default).
This makes sure that for each time step the accelerations from the previous time step are used.
This can be done via the GUI or using Python commands available in Abaqus similar to the following:

```python
from step import *
step1 = my_model.ImplicitDynamicsStep(name='Step-1', previous='Initial', timePeriod=1, nlgeom=ON, maxNumInc=1, haftol=1, initialInc=1, minInc=1, maxInc=1, amplitude=RAMP, noStop=OFF, nohaf=ON, initialConditions=OFF, timeIncrementationMethod=FIXED, application=QUASI_STATIC)
step1.Restart(frequency = 99999, overlay = ON)
```

The second command enables writing of restart files by Abaqus, which is required for running unsteady cases. This can be done from the GUI when the "Step" module is active in the viewport, by selecting "Output" in the top menu and subsequently "Restart Requests". *Frequency* should be put on *99999*, *overlay* *activated* (this spares disk space since only the last increment is kept) and _interval_ on _1_ (also see [this Abaqus documentation page](https://abaqus-docs.mit.edu/2017/English/SIMACAECAERefMap/simacae-t-simodbrestart.htm)).
Note that the step type "ImplicitDynamicsStep" is an example, it depends on the analysis procedure chosen to run the Abaqus part of the FSI simulation. For a steady simulation this could for example be "StaticStep":

```python
step1 = my_model.StaticStep(name='Step-1', previous='Initial', timePeriod=1.0, initialInc=1, minInc=1e-4, maxNumInc=10, nlgeom=ON, amplitude=RAMP)
```

The Abaqus wrapper tries to check if the increments comply with the `delta_t` setting and if needed adjusts it accordingly (notifying the user by raising a warning). The lines in the base-file (.inp) can look similar to this:

```
*Step, name=Step-1, nlgeom=YES, inc=1
*Dynamic,application=QUASI-STATIC,direct,nohaf,initial=NO
0.0001,0.0001,
```

It is required to have `initial=NO` when using the application transient fidelity (default when no application specified). This corresponds to the `initialConditons=OFF` setting when creating a step using Python. As mentioned earlier, this makes sure that for each time step the accelerations from the previous time step are used. The time step (0.0001) will in this case be replaced by settings found in the JSON file. More information for dynamic cases can be found on [this Abaqus documentation page](https://abaqus-docs.mit.edu/2017/English/SIMACAEKEYRefMap/simakey-r-dynamic.htm), for static cases on [this page](https://abaqus-docs.mit.edu/2017/English/SIMACAEKEYRefMap/simakey-r-static.htm).

### Note about choosing `ModelParts`
The created "surfaces" and "node sets" for load input and displacement output respectively, correspond to `ModelParts` in the CoCoNuT code, a representation of the data used for the coupling. It is strongly advised to sub-divide to fluid-structure interaction interface intelligently, depending on the geometry. As a rule of thumb it can be said that a surfaces at two sides of a sharp corner should be assigned to a different `ModelPart`. As the interpolation is based on shortest distance, issues can arise at sharp corners. Those are avoided by having different `ModelParts` at each side of the corner. Another reason to do this is because the code cannot handle elements with two or more faces being part of the same `ModelPart`. This situation would occur if the surface contains corners. An example is an airfoil where the suction side and pressure side belong to the same `ModelPart`: elements at the trailing edge will have (a) face(s) at both the pressure side and suction side. Even when the code would allow this, interpolation mistakes become likely, as a geometrical node or load point on the suction side could have a nearest neighbour on the pressure side, causing that the wrong data is used for interpolation.

## Log files

A general event log of the procedure can be found in the working directory, in a file named *`abaqus.log`*. For more detailed information on a certain time step, the .msg file written by Abaqus can be consulted. In CoCoNuT these are structured as follows: *`CSM_TimeA.msg`*, *`A`* being the time step. Typically multiple coupling iterations are done within each time step, so these .msg-files get overwritten by each new coupling iteration in the same time step.

## Restarting a calculation

For a restart of a calculation, say at time step `B`, it is necessary that all the Abaqus [simulation files](#files-written-in-the-working-directory-during-solve_solution_step) of the previous calculation at time step `B` are still present.
This is in accordance with the `save_restart` parameter defined in a higher [component](../../coupled_solvers/coupled_solvers.md). 
Particulary for the Abaqus solver wrapper, it is important that the files *`CSM_Time0Cpu0SurfaceAFaces.dat`*, *`CSM_Time0SurfaceAElements.dat`* and *`CSM_Time0SurfaceANodes.dat`* are still present from previous calculation.
The files *`CSM_Time0.inp`* and *`CSM_Restart.inp`* will be generated during initialization of the restarted calculation. 
This allows the user to alter some parameters in the input file before restart, e.g. altering output requests, boundary conditions or applying additional loads. 
Since Abaqus only uses the *`CSM_Restart.inp`* (which does not contain any mesh information) and output files of the previous calculation, it is pointless to change the mesh before a restart.

## End of the calculation

Unlike for [other solvers](../../coupling_components.md), CoCoNuT does not keep track of the time it takes to write case files.
The reason for this is that Abaqus starts and stops every single coupling iteration, see [here](#the-solve_solution_step-method) and therefore needs data files to be saved every time step in order to continue the calculation in a next coupling iteration or time step. 
Saving data files is as such seen as inherent to solving the solution step.

## Solver coupling convergence

The convergence criterion [solver coupling convergence](../../convergence_criteria/convergence_criteria.md#solver-coupling-convergence) can not be used for the Abaqus solver wrapper.

## Version specific documentation

### v2023
No major changes.

### v2024
Abaqus is now using Python 3.10 instead of Python 2.7.
