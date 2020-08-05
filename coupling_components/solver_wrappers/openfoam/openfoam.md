# SolverWrapperOpenFOAM

This is the documentation for all OpenFOAM solver-wrappers.
Currently only FSI simulations are supported, no other multiphysics problems.


## Parameters

This section describes the parameters in the JSON file, listed in alphabetical order.

parameter|type|description
---:|:---:|---
`case_file`|string|Name of the case file. It must be present in the folder `working_directory`.
`cores`|int|Number of processor cores to use when running Fluent (tested only on single node so far).
`dimensions`|int|Dimension used in flow solver: 2 for 2D and axisymmetric, 3 for 3D. 
`delta_t`|double|Fixed timestep size in flow solver. This parameter is usually specified in a higher `Component`.
`flow_iterations`|int|Number of non-linear iterations in Fluent per coupling iteration.
`fluent_gui`|bool|Set to `true` to run Fluent with the graphical interface.
<nobr>`hybrid_initialization`</nobr>|bool|If `true`, the hybrid initialization is used in Fluent before the first timestep. If `false`, the standard initialization is used, which requires that adequate reference values have been set in advance by the user in the `case_file`.
`interface_input`|dict|Keys are names of `ModelParts` for Fluent nodes. Each name must be the concatenation of an entry from `thread_names` and "_nodes". The values are (lists of) names of `Variables`.
`interface_output`|dict|Analogous to `interface_input`, but for Fluent faces ('_faces').
`max_nodes_per_face`|int|This value is used to construct unique IDs for faces, based on unique IDs of nodes (provided by Fluent). It should be at least as high as the maximum number of nodes on a face on the interface. Use e.g. 4 for rectangular faces, 3 for triangular faces.
`save_iterations`|int|Number of timesteps between consecutive saves of the Fluent case and data files.
`thread_names`|list|List with Fluent names of the surface threads on the interface. IS THE ORDER IMPORTANT? TODO
`timestep_start`|int|Timestep number to (re)start a transient FSI calculation. If 0 is given, the simulation starts from the `case_file`, else the code looks for the relevant case and data files. This parameter is usually specified in a higher `Component`.
`unsteady`|bool|`true` for transient FSI, `false` for steady FSI.  
`working_directory`|string|Absolute path to the working directory or relative path w.r.t the current directory.


`timestep_start` and `delta_t` are necessary parameters, but are usually defined in a higher `Component`. However, they can also be given directly as parameter of the solver-wrapper (e.g. for standalone testing). If they are defined both in a higher `Component` and in the solver-wrapper, then the former value is used and a warning is printed.

If different parameters are used with different Fluent versions, this should be specified both in this section and in the version specific documentation section.


## Overview of operation

The solver-wrapper itself consists of only one Python-file: `OpenFOAM_X.py`, with `X`the identifier of the OpenFOAM-version (e.g. `41` for OpenFOAM 4.1). Aside from this, the Python-file constructs a number of files such as `controlDict`, `pointDisplacement` and `decomposeParDict` by starting from the `_raw`-files in the solver-wrapper directory and replacing the names of the settings with their user-defined values (listed in the table above). Finally, using an OpenFOAM-solver in CoCoNuT requires the adaptation of the solver to accomodate the internal messaging system used during the FSI-simulation. Currently, only `pimpleFoam` and `interFoam` have been adapted; the solvers called by the solver-wrapper have the same name as the original OpenFOAM-solver but with `CoCoNuT_` added to the name. Upon using the OpenFOAM-wrapper, the solver to be used upon execution needs to be compiled using the default OpenFOAM-compilation method (loading the OpenFOAM-module and using `wmake`).

### The `__init__` method

During initialization, the journal and UDF files are adapted (parameter values are filled in) and copied to the `working_directory`. Fluent is then started in that directory using the parameters `cores`, `dimensions` and `fluent_gui`. Fluent then writes a case summery, so that `SolverWrapperFluentX` can link the interface thread names (specified in `thread_names`) to Fluent thread IDs, for use in the UDFs. Then the `Model` and `ModelParts` are created, based on data written by Fluent in timestep 0. After a restart, this same data must be found, i.e. if that file is removed from the folder `working_directory`, the simulation cannot be restarted. If the simulation was restarted, the coordinates `X`, `Y`, `Z` in the `Nodes` are updated to the current values. Finally, the `Interfaces` are made.

### Files created during simulation

In the file conventions, `A` is the timestep number and `B` the Fluent thread ID.

-   Fluent case and data files are saved as files of the form `case_timestepA.cas` and `case_timestepA.dat`.
-   Current node and face coordinates are passed from Fluent to Python with files of the form `nodes_timestepA_threadB.dat` and `faces_timestepA_threadB.dat`. 
-   The new node coordinates are passed from Python to Fluent with files of the form `nodes_update_timestepA_threadB`.
-   Pressure and tractions are passed from Fluent to Python with files of the form `pressure_traction_timestepA_threadB.dat`.
-   Files with extension `.coco` are used to pass messages between Python and Fluent (via the journal). 



## Setting up a new case

Following items should be set up and saved in the `case_file` (this list may be non-exhaustive):

-   additional UDFs must be configured
-   steady/unsteady (should match with the `unsteady` parameter)
-   2D, 3D or axisymmetric (should match with the `dimensions` parameter)
-   dynamic mesh for all zones, except the FSI interfaces
-   boundary conditions, material properties, numerical models, discretization schemes, operating conditions, turbulence modeling, convergence criteria
-   if `hybrid_initialization` is `false`, then defaults should be set for standard initialization

Following items are taken care of by CoCoNuT, and must therefore not be included in the `case_file`:

-   dynamic mesh for the FSI interfaces (which are defined in `thread_names`)
-   the timestep (`delta_t`) 
-   initialization of the solution field



## Version specific documentation

### 2019R1 (19.3)

This is currently the only version, so this section is still empty.

