# OpenFOAM

This is the documentation for all OpenFOAM solver wrappers. Currently, only FSI simulations are supported, no other
multiphysics problems.

## Parameters

This section describes the parameters in the JSON file, listed in alphabetical order.

parameter|type|description
---:|:---:|---
<nobr>`application`</nobr>|str|Name of the (adapted) OpenFOAM-solver to be used for the flow problem. This name should start with `coconut_`.
<nobr>`boundary_names`</nobr>|list| List of names of the patches corresponding to the interface. These names should match the patch names defined in the OpenFOAM-case.
`compile_clean`|bool|(optional) Default: `False`. If set to True, the adapted application will first clean and then compile.
`debug`|bool|(optional) Default: `False`. For every iteration, additional files are saved containing information on the input and output data of the solver.
`delta_t`|double|(optional) Fixed timestep size in flow solver.
`density`|int|The applied density of the fluid in an incompressible case. The density will be applied if `is_incompressible` is set to `true`.
<nobr>`interface_input`</nobr>|dict| List of dictionaries that describes the input `Interface`. This provides the  interface boundary conditions for the OpenFOAM solver. Each entry in the list has two keys: `model_part` and `variables`, with values as name of the model part and list of input variables, respectively. The input variables in the list should be chosen from the  `variables_dimensions` `dict` in  the file *`coconut/data_structure/variables.py`*. The model part name must be the concatenation of an entry from `boundary_names` and the string `_input`.
<nobr>`interface_output`</nobr>|dict| Analogous to `interface_input`, but here the name must be the concatenation of an entry from `boundary_names` and the string `_output`. The entries in the list provides boundary conditions for the other solver(s) participating in the coupled simulation.
<nobr>`is_incompressible`</nobr>|bool| Set it to `True` if the solver solves incompressible Navier-Stokes equation.
<nobr>`parallel`</nobr>|bool| Set it to `True` if OpenFOAM solver is required to run in parallel. The required decomposition method and number of cores should be provided in the *`<case_directory>/system/decomposeParDict`* file.
<nobr>`residual_variables`|list|(optional) A list containing OpenFOAM variables whose residuals you need to output. If provided, this will output the last initial residual of the pimple iterations for each FSI-coupling iteration in *`<case_directory>/residuals.csv`*.
<nobr>`time_precision`</nobr>|int|Number of digits after the decimal sign to be used in the name of the time step directories which are made during execution.
`timestep_start`|int|(optional) Time step to (re)start a transient FSI calculation from. If 0 is given, the simulation starts from t = 0, else the code looks for the relevant case and data files.
<nobr>`working_directory`</nobr>|str|Directory where the OpenFOAM-case is defined (and which contains the JSON-file).

`timestep_start` and `delta_t` are necessary parameters, but are usually defined already in the parameters of the
coupled solver. However, they can also be given directly as parameter of the solver wrapper (e.g. for standalone
testing). If they are defined both in the coupled solver and in the solver wrapper, then the former value is used, and a
warning is printed.

## Overview of the OpenFOAM-wrapper

The solver wrapper can be found in *`solver_wrappers/openfoam`* and consists of a python-file named *`vX.py`*, with `X`
the identifier of the OpenFOAM-version (e.g. *`v41.py`* for OpenFOAM 4.1). Finally, using an OpenFOAM-solver in CoCoNuT
requires the adaptation of the solver to accommodate communications with the `coupled_solvers` during the
FSI-simulation. Currently, only `pimpleFoam` and `interFoam` have been adapted; the solvers called by the solver wrapper
have the same name as the original OpenFOAM-solver but with the prefix `coconut_`. Upon using the OpenFOAM-wrapper, the
solver to be used upon execution needs to be compiled using the default OpenFOAM-compilation method (loading the
OpenFOAM-module and using `wmake` in the directory *`solver_wrapper/coconut_<solver_name>`*).

### The `__init__` method

During initialization, the case settings as defined in the JSON-file are read into the corresponding python-variables.
Afterwards, the required OpenFOAM-files are modified, and some checks are in place to verify whether the correct modules
are loaded. Next, the model-parts for the `interface_input/output` are created. Finally, the variables on each interface
are stored in the data-structure of CoCoNuT.

### The `initialize` method

The OpenFOAM-solver, which should be adapted to work with CoCoNuT, is launched in a subprocess. Other variables, such as
time step index and physical time, are initialized.

### The `initialize_solution_step` method

This function is called at the start of every time step. The time step index is increased by one, the iteration index is
set to zero, and the OpenFOAM file directory is made in which the time step data needs to be stored.

### The `solve_solution_step` method

The interface displacement is converted into an OpenFOAM-readable format (with the function `write_node_input`), after
which the OpenFOAM-loop is called in a subprocess. After completion of the OpenFOAM-loop, the
function `read_node_output` is called, which reads the interface loads from the corresponding OpenFOAM-directory (in
the *`postProcessing`* folder).

### The `finalize_solution_step` method

In this method, it is checked whether the other flow variables need to be written in an OpenFOAM-directory. This
save-option is done in OpenFOAM.

### The `finalize` method

The OpenFOAM-subprocess, which was launched in the `Initialize` method, is killed.

## Comments

- Files with extension `.coco` are used to pass messages between python-script and OpenFOAM. After sending a message
  from python to OpenFOAM, the python-script is paused until it receives the corresponding message from OpenFOAM. The
  OpenFOAM-solver is operating in an infinite `while`-loop until it receives a message from Python. Upon receipt of this
  message, OpenFOAM executes the corresponding action after which it sends a response to Python and reverts into its
  infinite loop (waiting mode).
- The aforementioned messaging procedure implies that OpenFOAM is constantly running during execution of the
  CoCoNuT-simulation. It is closed by the Python-code only in the `Finalize` method, but if CoCoNuT crashes, the
  OpenFOAM-programme keeps running. The user should take care to kill that OpenFOAM-loop manually (using `kill`
  or `pkill` in the Linux-terminal, e.g. `pkill coconut_*`).
- The interface displacement is stored in a `pointDisplacement` field, which is read in by OpenFOAM in every iteration (
  this required some adaptation of the solver, see next section). The dynamic mesh motion is handled by OpenFOAM itself.
- The interface loads are stored in the directory *`postProcessing`*, under the names *`PRESSURE_<boundary_name>`*
  and *`TRACTION_<boundary_name>`*. These are constructed from the file *`controlDict`*, defined in the `_init_` method
  of the solver wrapper in python.

## Overview of an OpenFOAM-solver used in CoCoNuT

Default OpenFOAM-solvers cannot be used in the CoCoNuT-framework, but need to be modified. The adapted solvers are
stored in the *`coconut/coupling_components/solver_wrappers`* directory and receive the name *`coconut_<solver_name>`*,
where `solver_name` is the name of the original solver, e.g. `pimpleFoam`. If a new solver needs to be adapted to work
with CoCoNuT, one of the already established solvers can work as an example. In brief, the following steps should be
undertaken:

- Except for the initial `include`-statements, the entire solver code should be put in a loop that starts
  with `while(true)`.
- Just before this while-loop, add the statement `runTime.run();`! This is important as it creates and initializes the
  functionObjects in the *`controlDict`* file which will be used for storing the interface loads.
- Due to the fact that the infinite loop is driven by `while(true)` and not the default
  OpenFOAM-expression `while(runTime.run())`, OpenFOAM does not check itself whether the final time step is reached. In
  proper CoCoNuT-operation, OpenFOAM is killed when the solver wrapper reaches the final time step.
- Inside the while-loop, a sleep-command is added such that the solver is not constantly checking the conditional
  statements.
- The while-loop contains several conditional statements, each of which check whether the python-code in CoCoNuT has
  sent a message to the OpenFOAM-solver. This message is sent by creating an empty file with a specific name in the
  OpenFOAM-directory. The following file names should be checked by the OpenFOAM-solver: *`next.coco`*
  , *`continue.coco`*, *`save.coco`*, *`stop.coco`*.
- If the file *`next.coco`* exists, the runTime-object should be increased by one. OpenFOAM should create a
  file *`next_ready.coco`* upon completion. Do not forget to delete the original *`next.coco`* file, which is advised to
  do just before creating the *`next_ready.coco`*, so near the end of the `if`-clause.
- If the file *`continue.coco`* exists, the flow equations need to be solved. This `if`-statement consequently contains
  most of the original solver definition, in which the flow equations are called in the same order as in the original
  CFD solver. OpenFOAM should create a file *`continue_ready.coco`* upon completion. Do not forget to delete the
  original `continue.coco`-file, which is advised to do just before creating the `continue_ready.coco`, so near the end
  of the `if`-clause.
- If the file *`save.coco`* exists, OpenFOAM checks whether the flow variables should be stored in corresponding files,
  according to the user-defined save interval. OpenFOAM should create a file *`save_ready.coco`* upon completion. Do not
  forget to delete the original *`save.coco`* file, which is advised to do just before creating the *`save_ready.coco`*,
  so near the end of the `if`-clause.
- If the file *`stop.coco`* exists, a `break`-statement should end the infinite loop (the subprocess is also killed in
  the python-code). OpenFOAM should create a file *`stop_ready.coco`* before breaking the `while`-loop. Do not forget to
  delete the original *`stop.coco`* file, which is advised to do just before creating the *`stop_ready.coco`*, so near
  the end of the `if`-clause.

## Setting up a new case

Following items should be present in the OpenFOAM-directory prior to launching CoCoNuT:

- The entire framework of the CFD-case in OpenFOAM which is to be used in the CoCoNuT simulation (so it should contain
  the `constant` and `system` directory as well as the *`0`* directory). The working directory should be defined as if
  you would like to run it as a CFD case. The working directory is defined in a JSON-file and therefore the
  CoCoNuT-files do not need to be in the same folder as the OpenFOAM-case.
- a JSON-file containing all the settings stipulated above.

Following files are modified or used as a reference to create files by CoCoNuT, and must be included in the original
OpenFOAM-directory:

- *`system/controlDict`* with compulsory arguments:

  key|value
      :---:|:---
  `writeControl`|`timeStep`
  `writeInterval`|required write interval
  `writeFormat`|`ascii`
  `timeFormat`|`fixed`
  `timePrecision`|required time precision based on `delta_t` (used in the names of time step folders)
  `runTimeModifiable`|`false`
  `adjustTimeStep`|`no`

- *`constant/dynamicMeshDict`*which contains the settings for OpenFOAM's dynamic motion solver
- *`system/decomposeParDict`* with the necessary decomposition of the fluid domain (if `cores`>1)
- *`0/pointDisplacement`* with all the boundary conditions, including `fixedValue` boundary condition for the FSI
  boundaries. This is used as a template for the *`pointDisplacementTmp`* to supply displacement boundary conditon (
  from structural solver) for the FSI-interface.

### Comments

- It is probably best to derive a new case from the directory containing FSI simulation with OpenFOAM
  in *`coconut/examples/`* in order to copy its structure.
- If you do not use an OpenFOAM-solver which is already converted for operation in CoCoNuT, you will need to adapt the
  solver yourself. This can be done in a rather straightforward way by taking a look at already implemented solvers.
  
## Version specific documentation

### v41 (OpenFOAM 4.1)

First version.

### v8 (OpenFOAM 8)

In this OpenFOAM version, the name of some commands and classes are different. The required changes are visible in the
version specific solver wrapper python file _`v8.py`_.

Moreover, there are some changes in the case setup. The following list provides some but is not all extensive:

- The interpolation point for `arc` definitions using `blockMesh` has to be on the arc and not just on its 'elongation'
- Using a macro expansion, such as `$R`, in combination with `#calc`, where an operator follows the macro expansion,
  requires using protection: use `#calc "${R}/2` instead of `#calc "$R/2`
- The file _`constant/turbulenceProperties`_ has been renamed to _`constant/momentumTransport`_
- The keyword `residualControl` in the `PIMPLE` dictionary is now called `outerCorrectorResidualControl`. The
  keyword `residualControl` still exists, but has a different meaning and is used as in `SIMPLE`.
- The file _`system/fvSolution`_ requires dictionaries `pcorr` and `pcorrFinal`, if the keyword `correctPhi` is `true`
  in _`system/fvSolution`_, which it is by default. The tolerance settings should be sufficiently low (similar to `p`).