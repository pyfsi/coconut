# Kratos- Structural

This is the documentation for the solver wrapper of the `StructuralMechanicsApplication` in Kratos.

## Parameters

This section describes the parameters in the JSON file, listed in alphabetical order.

parameter|type|description
---:|:---:|---
`cores`|int|Number of processor cores to use when running Kratos (works with 1 core, multi-processing is work in progress).
`delta_t`|double (optional)|Fixed timestep size in structural solver.
`input_file`|string| Project parameters file used by Kratos in json format.
`interface_input`|dict| List of dictionaries that describes the input `Interface`. This provides the  interface boundary conditions for the Kratos solver. Each entry in the list has two keys: `model_part` and `variables`, with values as name of the model part and list of input variables, respectively. The input variables in the list should be chosen from the  `variables_dimensions` `dict` in  the file *`data_structure/variables.py`*. The model part name must be the concatenation of an entry from `kratos_interface_sub_model_parts_list` and the `string` "_input".
`interface_output`|dict|Analogous to `interface_input`, but here the name must be the concatenation of an entry from `kratos_interface_sub_model_parts_list` and the `string` "_output". The entries in the list provides boundary conditions for the other solver(s) participating in the coupled simulation.
`solver_load_cmd`|string| Bash commmands for loading required modules and environmental variables to run Kratos,
`kratos_interface_sub_model_parts_list`|string| Names of sub-model parts used for input and output in Kratos,
`timestep_start`|int (optional)|Timestep number to (re)start a transient FSI calculation. If 0 is given, the simulation starts from t = 0, else the code looks for the relevant case and data files.  
<nobr>`working_directory`</nobr>|string|Absolute path to the working directory or relative path w.r.t the current directory.


`timestep_start` and `delta_t` are usually defined already in the parameters of the `coupled_solver`. However, they can also be given directly as a parameter of the solver-wrapper (e.g. for standalone testing). If they are defined both in the coupled solver and in the solver-wrapper, then the former value is used and a warning is printed.

If different parameters are used with different Kratos versions, this should be specified both in this section and in the version specific documentation section.


## Overview of operation

The solver-wrapper consists of 2 files, where `X` is the Kratos version without decimal, e.g. for version `6.0` `60`:

-   *`vX.py`*: defines the `SolverWrapperKratosStructureX` class
-   *`run_kratos_structural_X.py`*: The python file which runs Kratos in the background. This interacts with CoCoNuT for coupling.

### The `__init__` method

During initialization, the *`ProjectParameters.json`* file required for Kratos simulation is read and modified (parameter values are filled in) in the `working_directory`. One additional parameter called  `interface_sub_model_parts_list` is added in the *`ProjectParameters.json`* that tells Kratos about the interface model parts used in the coupling.  Kratos structural simulation is then started in that directory using the parameter `cores` (multi-processing not implemented yet). The solver wrapper waits for the Kratos simulation to output interface sub-model parts nodes, so that `SolverWrapperKratosStructureX` can create model parts in CoCoNuT for each entry in `interface_input` and `interface_output`. Finally, the interfaces are created.

### Files created during simulation


-   The interface sub-model parts nodes are saved as *`<sub_model_part_name>_nodes.csv`*.
-   The displacement from Kratos is written in a file named *`<sub_model_part_name>_displacement.csv`*.
-   Pressure and tractions are passed from python to Kratos with files of the form *`<sub_model_part_name>_pressure.csv`* and *`<sub_model_part_name>_surface_load.csv`*, respectively.
-   Files with extension ".coco" are used to pass messages between Python and Kratos. 



## Setting up a new case

Following items should be set up and saved in the `working_directory` (this list may be non-exhaustive):

-   *`ProjectParameters.json`* with all the required parameters and `end_time` set to a very high value, e.g. 1e9.
-   Mesh file with extension ".mdpa".
-   *`StructuralMaterials.json`* with the material properties.

Following items are taken care of by CoCoNuT, and must therefore will be automatically changed at the begining of the simulation:

-   the start_time (`timestep_start`)
-   the timestep (`delta_t`) 
-   initialization of the solution field



## Version specific documentation

### v60 (6.0)

First version.

