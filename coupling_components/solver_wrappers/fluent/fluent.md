# Fluent

This is the documentation for all Fluent solver wrappers.
The focus of this page is on fluid-structure interaction (FSI). Currently no other multiphysics problems are supported, but these are envisioned in future developments. FSI with inviscid flows is not supported, but that would be straightforward to add if required.


## Fluid-structure interaction with Fluent

ANSYS Fluent can be used to solve the flow field in partitioned FSI simulations. The FSI interface consists of one or several *surface threads* in Fluent. 
During each FSI iteration, the solver wrapper imposes a certain displacement on this FSI interface. Fluent then automatically deforms the rest of the mesh (using dynamic mesh capabilities), solves the flow field and returns the calculated loads on the FSI interface to the solver wrapper. 
Displacements are applied in Fluent in *nodes*. As displacements are the input of Fluent, these nodes are collected in one or more `ModelParts` in the *input* `Interface`. 
Loads (pressure and traction) can be obtained from the Fluent *faces*. As loads are the output of Fluent, these faces are collected in one or more `ModelParts` in the *output* `Interface`.
A fixed naming convention is used for the Fluent `ModelParts` in CoCoNuT: each `ModelPart` gets the name of the corresponding Fluent surface thread, plus "_nodes" or "_faces", depending on their content. As a consequence, each `ModelPart` in the input `Interface` has a counterpart in the output `Interface`. 
More information about `ModelParts` and `Interfaces` can be found in the [data structure documentation](../../../data_structure/data_structure.md).


## Parameters

This section describes the parameters in the JSON file, listed in alphabetical order.

parameter|type|description
---:|:---:|---
`case_file`|string|Name of the case file. It must be present in the folder specified by `working_directory`. The corresponding data file must also be present but has no key in the JSON file.
`cores`|int|Number of processor cores to use when running Fluent (tested only on single node so far).
`dimensions`|int|Dimension used in flow solver: `2` for 2D and axisymmetric, `3` for 3D. 
`delta_t`|float|Fixed time step size in flow solver. This parameter is usually specified in a higher `Component`.
`flow_iterations`|int|Number of non-linear iterations in Fluent per coupling iteration.
`fluent_gui`|bool|Set to `true` to run Fluent with the graphical interface.
`interface_input`|list|List of dictionaries to describe the input `Interface` (Fluent nodes). Each dictionary defines one `ModelPart` with two keys: `model_part` contains the name of the `ModelPart` and `variables` contains a list of variable names. Each `ModelPart` name must be the concatenation of an entry from `thread_names` and "_nodes". The variable names must be chosen from *`data_structure/variables.py`*. 
`interface_output`|dict|Analogous to `interface_input`, but for the output `Interface` (Fluent faces). Each `ModelPart` name must be the concatenation of an entry from the file `thread_names` and "_faces".
<nobr>`max_nodes_per_face`</nobr>|int|This value is used to construct unique IDs for faces, based on unique IDs of nodes (provided by Fluent). It should be at least as high as the maximum number of nodes on a face on the interface. Use e.g. 4 for rectangular faces, 3 for triangular faces.
`save_iterations`|int|Number of time steps between consecutive saves of the Fluent case and data files.
`thread_names`|list|List with Fluent names of the surface threads on the FSI interface. 
`timestep_start`|int|Time step number to (re)start a transient FSI calculation. If 0 is given, the simulation starts from the `case_file`, else the code looks for the relevant case and data files. This parameter is usually specified in a higher `Component`.
`unsteady`|bool|`true` for transient FSI, `false` for steady FSI.
`working_directory`|string|Absolute path to the working directory or relative path w.r.t the current directory.


`timestep_start` and `delta_t` are necessary parameters, but are usually defined already in the parameters of the coupled solver. However, they can also be given directly as parameter of the solver wrapper (e.g. for standalone testing). If they are defined both in the coupled solver and in the solver wrapper, then the former value is used and a warning is printed.

If different parameters are used with different Fluent versions, this should be specified both in this section and in the version specific documentation section.


## Overview of operation

The solver wrapper consists of 3 files (with X the Fluent version, e.g. "v2019R1"):

-   *`X.py`*: defines the `SolverWrapperFluentX` class, 
-   *`X.jou`*: Fluent journal file to interactively run the FSI simulation, written in Scheme, 
-   *`X.c`*: Fluent UDF file that implements additional functionality used in Fluent, written in C.

### The `__init__` method

During initialization, the journal and UDF files are adapted (parameter values are filled in) and copied to the `working_directory`. Fluent is then started in that directory using the parameters `cores`, `dimensions` and `fluent_gui`. Fluent then writes a case summary, so that the solver wrapper can link the interface thread names (specified in `thread_names`) to Fluent thread IDs, for use in the UDFs. Then the `Model` and `ModelParts` are created, based on data written by Fluent in time step 0. After a restart, this same data must be found, i.e. if those data files are removed from the `working_directory` folder, the simulation cannot be restarted. Finally, the `Interfaces` are created. 

### Files created during simulation

In these file conventions, A is the time step number and B the Fluent thread ID.

-   Fluent case and data files are saved as files of the form *`case_timestepA.cas`* and *`case_timestepA.dat`*.
-   Current node and face coordinates are passed from Fluent to CoCoNuT with files of the form *`nodes_timestepA_threadB.dat`* and *`faces_timestepA_threadB.dat`*. 
-   The new node coordinates are passed from CoCoNuT to Fluent with files of the form *`nodes_update_timestepA_threadB.dat`*.
-   Pressure and tractions are passed from Fluent to CoCoNuT with files of the form *`pressure_traction_timestepA_threadB.dat`*.
-   Files with extension *`.coco`* are used to exchange messages between CoCoNuT and Fluent. 



## Setting up a new case

Following items should be set up and saved in the Fluent case file (this list may be non-exhaustive):

-   additional UDFs must be configured, 
-   steady/unsteady (should match with the `unsteady` parameter),
-   2D, 3D or axisymmetric (should match with the `dimensions` parameter),
-   dynamic mesh zones for all deforming surfaces, except for the FSI interfaces,
-   dynamic mesh smoothing parameters,
-   boundary conditions, material properties, numerical models, discretization schemes, operating conditions, turbulence modeling, convergence criteria.

A data file should also be present with the fields either initialized or containing the results of a previous calculation.

Following items are taken care of by CoCoNuT, and must therefore not be included in the Fluent case file:

-   dynamic mesh zones for the FSI interfaces (these are defined in `thread_names`),
-   the time step (`delta_t`).


## Version specific documentation

### v2019R1 (19.3.0)

First version.

### v2019R2 (19.4.0)

No changes.

### v2019R3 (19.5.0)

The solutions in this version are (slightly) different because the *Rhie-Chow face flux interpolation in the pressure-based solver* has changed. This setting can be reverted with the TUI command `solve set previous undo-2019r3 y n`, which is included in *`v2019R3.jou`*.

The results can be slightly different when restarts are used for multi-core simulations for the following reason: *For parallel cases with smoothing that do not use dynamic load balancing, a zonal partitioning with Laplace smoothing will automatically be applied when the file is read, which should result in better load balancing for the mesh smoothing calculations.*
After a restart, the partitioning can be different and hence the mesh deformation can be slightly different. 

### v2020R1 (20.1.0)

Same behavior as v2019R3.