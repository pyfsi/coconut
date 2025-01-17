# AbaqusCSE

This is the documentation for the new version of the Abaqus solver wrapper, which uses the Abaqus Co-Simulation Engine (CSE) avoiding the need to start and stop Abaqus every coupling iteration.
Refer to [Abaqus](../abaqus.md) for the old version of the Abaqus solver wrapper.

Abaqus is a structural solver implementing the finite element method.
Currently, this wrapper only supports FSI simulations, no other multi-physics problems.

!!! failure "For two-dimensional cases, the traction forces are not taken into account"

!!! failure "Restart has not yet been implemented"

??? info "Terminology"

    - A _step_ in Abaqus is a convenient period of time (during which a load is applied for example). Within the context of CoCoNuT it usually coincides with the complete duration of the calculation. It should not be confused with the _time step_.
    - An _increment_ in Abaqus is usually equal to the time step of Abaqus, unless subcycling is applied within the Abaqus solver. Then a time step can be subdivided in multiple increments.
    - In an implicit or strongly coupled FSI calculation, the structural solver is called multiple times per time step. A call of the the flow solver and subsequent call of the structural solver is called a _coupling iteration_. All coupling iterations are part of the same attempt of the same increment (unless subcycling is used).

??? warning "Subcycling"

    Subcycling has not been tested for this solver wrapper. Note that the maximum number of increments is now set to the number of time steps. In order to apply subcycling, this needs to be modified.

??? success "Licensing"

    If the **Abaqus license server** needs to be specified explicitly, it is advised to do this in the [solver modules file](../../../README.md#checking-the-solver-modules) by setting the variable _LM_LICENSE_FILE_.
    The file _`abaqus_v6.env`_ is the environment file for the Abaqus solver and will be automically created by the CoCoNuT Python wrapper.

## Fluid-structure interaction with Abaqus
Abaqus (Dassault Systèmes) can be used to solve for the structural displacement/deformation in partitioned FSI-simulations.
The FSI interface consist of a _surfaces_ in the Abaqus model, where pressure and traction loads are applied.
The loads are applied in the element centers, the displacements are exported in the element nodes.
The input loads are collected in one or more `ModelParts` in the _input_ `Interface`,
the output nodes are collected in one or more `ModelParts` of the _output_ `Interface`.
Each `ModelPart` on the input `Interface` has a counterpart on the output `Interface`.
More information about `ModelParts` and `Interface` can be found in the [data structure documentation](../../../data_structure/data_structure.md).

## Parameters
This section describes the parameter settings in the JSON file.
It can be useful to have a look at a JSON file of one of the examples in the _`examples`_ folder.

|                                         parameter | type | description                                                                                                                                                                                                                                                                                                                                                                                                                |
|--------------------------------------------------:|:----:|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|                                           `cores` | int  | Number of cores to be used by Abaqus.                                                                                                                                                                                                                                                                                                                                                                                      |
|                                           `debug` | bool | (optional) Default: `false`. Additional files are created or retained that are helpful for debugging, such as Abaqus send-and-receive file (Abaqus.SRE) or the input and output data of the solver for every coupling iterations.                                                                                                                                                                                          |
|                                      `dimensions` | int  | Dimensionality of the problem (2 or 3).                                                                                                                                                                                                                                                                                                                                                                                    |
| <nobr>`disable_modification_of_input_file`</nobr> | bool | (optional) Default: `false`. If `true`, the input file is not modified by CoCoNuT but simply copied to a file named Abaqus.inp. The user is responsible to correctly define the co-simulation settings, the time step, the duration, the output and the restart behaviour.                                                                                                                                                 |
|                                 `interface_input` | list | Should contain a dictionary for each input `ModelPart` to be created, having a key `model_part` and a key `variables`. The order of the list should correspond to the order of the Abaqus _surfaces_ provided in the parameter `surfaces`, as well as to the order of the model parts in the other solver. The value of the key `variables` should be a list of the variables `pressure` and `traction`. The `interface_input` corresponds to the face centers of the elements on the surface.                  |
|                                `interface_output` | list | Similar to `interface_input` but contains the output `ModelParts` for Abaqus geometrical nodes. In this case the `variables` key specifies the output variable, chosen from _`data_structure/variables.py`_. Here, the value of the key `variables` should be a `displacement`. Again, the order should correspond to the `surfaces` as well as to the `Interface` definitions of the solver to which Abaqus is coupled. |
|                                      `input_file` | str  | Name of the Abaqus input file (.inp) provided by the user. <br> Example: _`case.inp`_                                                                                                                                                                                                                                                                                                                                      |
|                                            `port` | int  | (optional) If not provided, the operating system will choose a free port. This port is used for the communication between the AbaqusWrapper executable and the CSE.                                                                                                                                                                                                                                                        |
|                                    `save_results` | int  | (optional) Default: `1`. This will modify the frequency of storing fields or history in the _`.odb`_ file. If no output definition is present in the input file, one is created for the _PRESELECT_ variables.                                                                                                                                                                                                             |
|                                        `surfaces` | list | List of the surface names defined in Abaqus. Note that the order has to be the same as those of the interface definitions. Take the Abaqus naming conventions into account: surfaces names are capitalized unless specified within double quotes in the input file, and the name is often prepended for example with ASSEMBLY_ (“Assembly.Part” is translated to “Assembly_Part”).                                         |
|                               `working_directory` | str  | Relative path to the directory in which Abaqus will be executed and where all structural information will be stored. This folder is typically called _`CSM`_, but any name is allowed. It should be created before execution and contain the input file.                                                                                                                                                                   |

The following parameters are usually defined in the top-level settings of the JSON file, but they can also be given directly as parameter of the solver wrapper (e.g. for standalone testing).
If they are defined in both locations, a warning is printed and the top-level value is used.

|                          parameter | type  | description                                                                                                                                              |
|-----------------------------------:|:-----:|----------------------------------------------------------------------------------------------------------------------------------------------------------|
|                          `delta_t` | float | Size of the time step in Abaqus.                                                                                                                         |
| <nobr>`number_of_timesteps`</nobr> |  int  | The amount of time steps to run the calculation. For a steady calculation, the value should be `1`.                                                      |
|                     `save_restart` |  int  | Indicates the time step interval at which files for restart have to be saved. A minus sign indicates only the files from the last interval are retained. |
|                   `timestep_start` |  int  | Time step to start from. Data should be available at this time step. For a new simulation this value will typically be `0`.                              |


## Overview of operation
The communication with Abaqus is achieved through the Co-Simulation Engine (CSE).
Therefore, CoCoNuT will launch 3 processes simultaneously, as shown below.
The dashed arrows represent communication at the start of the calculations, while the full arrows represent communication of pressure, traction and displacement on run time.

![overview](images/abaqus_cse_structure.svg "Overview of operation of AbaqusCSE")

On the far right, the **Abaqus** process is shown, which is the actual solver.
The used mesh information and settings are contained within an input file _`Abaqus.inp`_.
CoCoNuT crates this file based on the provided `input_file` in the JSON file.
The time step, duration, number of increments, output frequency, restart frequency and co-simulation settings are all updated or added by the Python solver wrapper (unless `disable_modification_of_input_file` is `true`). This includes the load definition on the `ModelParts` and a basic output request.

The communication with this Abaqus process is taken care of by the **Co-Simulation Engine (CSE)**.
This is a second process launched by the Python solver wrapper.
Its settings are dictated by the file _`CSE_config.xml`_ which is also created by Python solver wrapper.

Finally, a third process is launched: the executable **AbaqusWrapper**.
Its C++ code is stored in _`coupling_components/solver_wrappers/abaqus_cse/AbaqusWrapper`_ and automatically compiled by the Python solver wrapper.
At the start of a calculation, it receives input written by the Python solver wrapper by reading the file _`AbaqusWrapper_input.txt`_.
When the calculation is running, the Python solver wrapper communicates with the AbaqusWrapper using text files such as: _`pressure_mp0.txt`_.
The suffix _mp0_ refers to model part 0, that is the first model part.

!!! success "Log files"

    Since the Python wrapper will launch three separate processes, it is wise to check the log file of each process should an error occur.

    - The Abaqus process is launched in the `working_directory`, and the log file is named _`abaqus.log`_. Usually, useful information can be found in _`Abaqus.msg`_ and occasionally _`Abaqus.dat`_.
    - The CSE process is launched in a subfolder, _`CSE`_, of the `working_directory`, and the log file is named _`cse.log`_.
    - The AbaqusWrapper process is launched in a subfolder, _`AbaqusWrapper`_,  of the `working_directory` , and the log file is named _`AbaqusWrapper.log`_.

    For debugging purposes, it is usefull to set the `debug` JSON parameter on true, resulting in the output of additional files.

## Setting up a case: creation of an Abaqus input file (.inp)
The Abaqus solver wrapper is configured to start from an input file which contains all necessary information for the calculation.
Its name should be specified in the JSON file via the parameter `input_file`, and it should be located in the `working_directory`.
For the remainder of this section this file will be referred to as "base file".

Creation of the base file is not considered a part of the solver wrapper functionality as it is case-specific.
However, in order for the solver wrapper to work, the base file has to comply with certain general conditions.
This section aims at informing the user about the requirements for the base file.

### General
The base file needs to be of the ".inp" type, this is an "input file for Abaqus". ".inp-files" are created via Abaqus by, after configuration, creating a "job" and requesting a "write input" for that job.
These files can be opened in Abaqus by using "file > import > model".
The base file has to contain all necessary information about the structural model, which includes:

 - Mesh defining the structure geometry and discretization.
    - Also, the element type needs to be defined.
 - Material properties.
 - Boundary conditions.
 - Surfaces where external loads need to be applied (one surface per `ModelPart`).
    - Here _"Surface"_ refers to nomenclature of the Abaqus software itself. The name can be found as such in the Abaqus working tree.
 - A _Step_ definition, which contains solver settings. The following type of analyses are some examples (it is advised to explicitly set them based on this documentation rather than leaving it to Abaqus to fill in a default):
    - Implicit dynamic, application quasi-static
    - Implicit dynamic, application moderate dissipation
    - Implicit dynamic, application transient fidelity
    - Static general
 - Additional loads not dependent on the flow solver.
 
Abaqus models contain parts and those parts are used to create assemblies.
The assembly, determines the position and orientation that will be used.

Abaqus has a GUI as well as a Python interface (which is also accessible) via the GUI. References to both the Python interface and GUI will be made below.

### Creating surfaces in Abaqus
Per surface in the fluid-structure interface (where loads and displacements need to be exchanged) a "surface" should be created. There are multiple possibilities to create these surfaces:

 - From the geometry: when the geometry has been defined in Abaqus itself, the geometry faces can easily be selected in the GUI. This method is often the most straightforward, but the Abaqus model should contain the geometry.
 - From the mesh: when the geometry is not available (this can for example be the case when a mesh has been imported), a surface can be defined by selecting multiple mesh faces. As a surface typically covers many mesh faces, it is useful there to select the regions "by angle", which uses the angle between mesh faces to determine whether adjacent faces should be selected. This way the surface selection can be extended until a sharp corner is met.
 - By converting a "node set" containing all the nodes on the surface and then calling the `SurfaceFromNodeSet` method which can be found in the _`make_surface.py`_ file in _`coupling_components/solver_wrappers/abaqus/extra`_. It is advised to copy this file to the working directory.

An example on the use of `SurfaceFromNodeSet` (via the Python console in Abaqus or a Python script for Abaqus):
```python
from make_surface import SurfaceFromNodeSet
my_model = mdb.models['Model-1']
my_assembly = my_model.rootAssembly  
my_instance = my_assembly.instances['PART-1-1']
inputSurfaceA = SurfaceFromNodeSet(my_assembly, my_instance, 'NODESET_NAME_A', 'SURFACE_NAME_A')
```

!!! tip "Note about choosing model parts and surfaces"

    The created "surfaces" for load input and displacement output, correspond to `ModelParts` in the CoCoNuT code. It is strongly advised to subdivide to fluid-structure interaction interface intelligently, depending on the geometry.
    As a rule of thumb it can be said that a surfaces at two sides of a sharp corner should be assigned to a different `ModelPart`.
    Issues can arise at sharp corners, as the interpolation is based on shortest distance.
    Those are avoided by having different `ModelParts` at each side of the corner.
    
    Another reason to do this is because the code cannot handle elements with two or more faces being part of the same `ModelPart`. 
    This situation would occur if the surface contains corners.
    An example is an airfoil where the suction side and pressure side belong to the same `ModelPart`: elements at the trailing edge will have (a) face(s) at both the pressure side and suction side.
    Even when the code would allow this, interpolation mistakes become likely, as a geometrical node or element center on the suction side could have (a) nearest neighbour(s) on the pressure side, leading to wrong data being used for interpolation.

## End of the calculation

Unlike for the [other solvers](../../coupling_components.md), CoCoNuT does not keep track of the time it takes to write output files, since this is controlled by Abaqus and triggered by CoCoNuT.
Nevertheless, the output frequency is controlled by CoCoNuT and can be set with the JSON file parameter `save_results`.

## Solver coupling convergence

The convergence criterion [solver coupling convergence](../../convergence_criteria/convergence_criteria.md#solver-coupling-convergence) has not yet been implemented for this Abaqus solver wrapper.

## Version specific documentation

### v2023
No major changes.

### v2024
Abaqus is now using Python 3.10 instead of Python 2.7.
