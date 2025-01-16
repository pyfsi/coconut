# KratosStructure

KratosMultiphysics is an open source framework for finite element simulations. More information on Kratos and the source code can be found on the [Kratos GitHub page](https://github.com/KratosMultiphysics). 
The solver wrapper of `StructuralMechanicsApplication` from Kratos has been implemented in CoCoNuT. This a documentation of this solver wrapper. 

## Parameters

This section describes the parameters in the JSON file, listed in alphabetical order.

|                                           parameter |  type  | description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
|----------------------------------------------------:|:------:|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|                                             `debug` |  bool  | (optional) Default: `false`. For every iteration, text files are saved with the input and output data of the solver.                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
|                                        `input_file` |  str   | Project parameters file used by Kratos in JSON format. In the [example cases](../../../examples/tube/fluent3d_kratos_structure3d/tube_fluent3d_kratos_structure3d.md), this is typically called *`ProjectParameters.json`*.                                                                                                                                                                                                                                                                                                                                                   |
|                                   `interface_input` |  dict  | List of dictionaries that describes the input `Interface`. This provides the  interface boundary conditions for the Kratos solver. Each entry in the list has two keys: `model_part` and `variables`, with values as name of the model part and list of input variables, respectively. The input variables in the list should be chosen from the  `variables_dimensions` `dict` in  the file *`coconut/data_structure/variables.py`*. The model part name must be the concatenation of an entry from `kratos_interface_sub_model_parts_list` and the string `_input`.         |
|                                  `interface_output` |  dict  | Analogous to `interface_input`, but here the name must be the concatenation of an entry from `kratos_interface_sub_model_parts_list` and the string `_output`. The entries in the list provides boundary conditions for the other solver(s) participating in the coupled simulation.                                                                                                                                                                                                                                                                                          |
| <nobr>`kratos_interface_sub_model_parts_list`<nobr> |  str   | Names of sub-model parts used for input and output in Kratos.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
|                               `pressure_directions` |  list  | A list containing 1 or -1 designating the direction in which the pressure is applied: a positive value means NEGATIVE_FACE_PRESSURE is used. This implies applying a pressure on the - face, which hence goes in the direction of the normal (-1,1), in other words 1 should be used when the normal points outwards from the fluid domain. A negative unit value results in the opposite directions, in other words -1 should be used when the normal points into the fluid domain. The length of the list must be equal to that of `kratos_interface_sub_model_parts_list`. |
|                        `print_coupling_convergence` |  bool  | (optional) Default `false`. If `true` and if the solver coupling convergence is checked, a statement is printed when the solver converges in the first solver iteration, see [solver coupling convergence](#solver-coupling-convergence).                                                                                                                                                                                                                                                                                                                                     |
|                                `residual_variables` |  list  | (optional) A list containing variables as reported in the log file (e.g. DISPLACEMENT, RESIDUAL or RESIDUAL DISPLACEMENT) whose residuals you need to output. If provided, this will output the last residual for each FSI-coupling iteration in *`<case_directory>/residuals.csv`*. For different element types, the names might be different and/or changes might be required to parse the log file. For the correct names see the Kratos log file.                                                                                                                         |
|                              `structure_iterations` |  int   | (optional) Maximum number of Newton iterations in Kratos per coupling iteration. If not provided, the value in `input_file` is used.                                                                                                                                                                                                                                                                                                                                                                                                                                          |
|                                    `timestep_start` |  int   | (optional)                                                                                                                                                                                                                                                                                                                                                                                               |
|                                 `working_directory` |  str   | Path to the working directory (i.e. where the `input_file` for Kratos is located), either absolute or relative w.r.t the current directory (i.e. from where the analysis is started).                                                                                                                                                                                                                                                                                                                                                                                         |


The following parameters are usually defined in the top-level settings of the JSON file, but they can also be given directly as parameter of the solver wrapper (e.g. for standalone testing).
If they are defined in both locations, a warning is printed and the top-level value is used.

|        parameter | type  | description                                                                                                                                                                         |
|-----------------:|:-----:|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|        `delta_t` | float | Fixed time step size in structural solver                                                                                                                                           |
|   `save_restart` |  int  | (optional) Default: `0`. Indicates the time step interval at which files for restart have to be saved. A minus sign indicates only the files from the last interval are retained.   |
| `timestep_start` |  int  | Time step to (re)start a transient FSI calculation from. If `0` is given, the simulation starts from t = 0, else the code looks for the relevan  files.                             |

If different parameters are used with different Kratos versions, this should be specified both in this section and in the version specific documentation section.


## Overview of operation

The solver wrapper consists of 2 files, where `X` is the Kratos version without decimal, e.g. for version 6.0 this becomes `60`:

-   *`vX.py`*: defines the `SolverWrapperKratosStructureX` class,
-   *`run_kratos_structural_X.py`*: The Python file which runs Kratos in the background. This interacts with CoCoNuT for coupling.

### The `__init__` method

During initialization, the *`ProjectParameters.json`* file required for Kratos simulation is read and modified (parameter values are filled in) in the `working_directory`. 
One additional parameter called  `interface_sub_model_parts_list` is added in the *`ProjectParameters.json`* that tells Kratos about the interface model parts used in the coupling. 
Kratos structural simulation is then started in that directory using the parameter `cores` (multiprocessing not implemented yet). 
The solver wrapper waits for the Kratos simulation to output interface sub-model parts nodes, so that `SolverWrapperKratosStructureX` can create model parts in CoCoNuT for each entry in `interface_input` and `interface_output`. 
Finally, the interfaces are created.

### Files created during simulation

-   The interface sub-model parts nodes are saved as *`<sub_model_part_name>_nodes.csv`*.
-   The displacement from Kratos is written in a file named *`<sub_model_part_name>_displacement.csv`*.
-   Pressure and tractions are passed from Python to Kratos with files of the form *`<sub_model_part_name>_pressure.csv`* and *`<sub_model_part_name>_surface_load.csv`*, respectively.
-   Files with extension *`.coco`* are used to pass messages between Python and Kratos.



## Setting up a new case

Following items should be set up and saved in the `working_directory` (this list may be non-exhaustive):

-   *`ProjectParameters.json`* with all the required parameters (see the [example cases](../../../examples/tube/fluent3d_kratos_structure3d/tube_fluent3d_kratos_structure3d.md))
-   Mesh file with extension *`.mdpa`*
-   *`StructuralMaterials.json`* with the material properties

Following items are taken care of by CoCoNuT, and must therefore will be automatically changed at the begining of the simulation:

-   The start time (`timestep_start`)
-   The time step (`delta_t`)
-   Initialization of the solution field

### Tip
````
When starting a simulation from time=0 after having restarted in the same folder, clean the directory or delete the case folder and copy the required files from a "set up" folder.
````
The reason for this is that the `ProjectParameters.json` file is modified when performing restart in such a way that it can no longer be used to start from time=0.

## Post-processing

Kratos `StructuralMechanicsApplication` v91 allows for saving VTK files. The number of files saved is set in `ProjectParameters.json` and is not modified by CoCoNuT.
Note that the zeroth instance corresponds to the first time step.


## Solver coupling convergence

The convergence criterion [solver coupling convergence](../../convergence_criteria/convergence_criteria.md#solver-coupling-convergence) is implemented for the Kratos StructuralMechanics solver wrapper.

When using the `residual_criterion` in Kratos as convergence criterion, it's important to set the relative tolerance to `-1`, effectively disabling it.
This is needed because Kratos sets the relative residual to zero in the first time step, wrongly leading to immediate convergence.
Although this is also the case when other coupling convergence criteria are used, it deserves special attention in case of the solver coupling convergence criterion.

Below, an example of the *`ProjectParameters.json`* file is included using the `and_criterion`.
This criterion combines the displacement and residual criteria and declares the solver converged only if both criteria are satisfied.
```json
{
  ...
  "convergence_criterion": "and_criterion",
  "displacement_relative_tolerance": 1e-7,
  "displacement_absolute_tolerance": 1e-7,
  "residual_relative_tolerance": -1,
  "residual_absolute_tolerance": 1e-7,
  ...
}
```

## Version specific documentation

### v94 (9.4)

Base version.

#### Kratos installation

This version supports installation using `pip`. To install only StructuralMechanicsApplication of Kratos, use the following command:
```bash
pip install KratosMultiphysics KratosStructuralMechanicsApplication KratosLinearSolversApplication
```
Note, the linear solvers in `KratosLinearSolversApplication` is sometimes used by the `StructuralMechanicsApplication`, therefore it is installed along with it.

To install all the available applications in Kratos, use the following command:
```bash
pip install KratosMultiphysics-all
```
