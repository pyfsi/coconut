# OpenFOAM

This is the documentation for all OpenFOAM solver wrappers.
Currently, the use as the flow solver in FSI simulations is supported, no other multiphysics problems.

## Parameters

This section describes the parameters in the JSON file, listed in alphabetical order.

|                                parameter |  type  | description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
|-----------------------------------------:|:------:|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|                            `application` |  str   | Name of the (adapted) OpenFOAM-solver to be used for the flow problem. This name should start with `coconut_`. For OpenFOAM 11, only `coconut_foamRun` is possible.                                                                                                                                                                                                                                                                                                                                                                              |
|                         `boundary_names` |  list  | List of names of the patches corresponding to the interface. These names should match the patch names defined in the OpenFOAM-case.                                                                                                                                                                                                                                                                                                                                                                                                              |
|                          `compile_clean` |  bool  | (optional) Default: `false`. If set to true, the adapted application will first clean and then compile.                                                                                                                                                                                                                                                                                                                                                                                                                                          |
|                                  `debug` |  bool  | (optional) Default: `false`. For every iteration, additional files are saved containing information on the input and output data of the solver.                                                                                                                                                                                                                                                                                                                                                                                                  |
|                                `delta_t` | double | Fixed timestep size in flow solver.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
|                                `density` | double | (optional) Density of the fluid in an incompressible case. The density is multiplied with the kinematic pressure and traction. Required if an incompressible application is used, such as coconut_pimpleFoam. More information can be found [here](#treatment-of-kinematic-pressure-and-traction).                                                                                                                                                                                                                                               |
|                        `interface_input` |  dict  | List of dictionaries that describes the input `Interface`. This provides the  interface boundary conditions for the OpenFOAM solver. Each entry in the list has two keys: `model_part` and `variables`, with values as name of the model part and list of input variables, respectively. The input variables in the list should be chosen from the  `variables_dimensions` `dict` in  the file *`coconut/data_structure/variables.py`*. The model part name must be the concatenation of an entry from `boundary_names` and the string `_input`. |
|                       `interface_output` |  dict  | Analogous to `interface_input`, but here the name must be the concatenation of an entry from `boundary_names` and the string `_output`. The entries in the list provides boundary conditions for the other solver(s) participating in the coupled simulation.                                                                                                                                                                                                                                                                                    |
|                               `parallel` |  bool  | Set it to `true` if OpenFOAM solver is required to run in parallel. The required decomposition method and number of cores should be provided in the *`<case_directory>/system/decomposeParDict`* file.                                                                                                                                                                                                                                                                                                                                           |
| <nobr>`print_coupling_convergence`<nobr> |  bool  | (optional) Default `false`. If `true` and if the solver coupling convergence is checked, a statement is printed when the solver converges in the first solver iteration, see [solver coupling convergence](#solver-coupling-convergence).                                                                                                                                                                                                                                                                                                        |
|                     `residual_variables` |  list  | (optional) A list containing OpenFOAM variables whose residuals you need to output. If provided, this will output the last initial residual of the pimple iterations for each FSI-coupling iteration in *`<case_directory>/residuals.csv`*.                                                                                                                                                                                                                                                                                                      |
|                         `time_precision` |  int   | Number of digits after the decimal sign to be used in the name of the time step directories which are made during execution.                                                                                                                                                                                                                                                                                                                                                                                                                     |
|                         `timestep_start` |  int   | Time step to (re)start a transient FSI calculation from. If 0 is given, the simulation starts from t = 0, else the code looks for the relevant case and data files.                                                                                                                                                                                                                                                                                                                                                                              |
|                      `working_directory` |  str   | Directory where the OpenFOAM-case is defined (and which contains the JSON-file).                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |

`timestep_start` and `delta_t` are necessary parameters, but are usually defined already in the parameters of the
coupled solver. However, they can also be given directly as parameter of the solver wrapper (e.g. for standalone
testing). If they are defined both in the coupled solver and in the solver wrapper, then the former value is used, and a
warning is printed.

## Setting up a new case

Following items should be present in the OpenFOAM-directory prior to launching CoCoNuT:

- The entire framework of the CFD-case in OpenFOAM which is to be used in the CoCoNuT simulation (so it should contain
  the `constant` and `system` directory as well as the *`0`* directory). The working directory should be defined as if
  you would like to run it as a CFD case. This `working directory` is defined in the JSON-file.
- The necessary parameters in the JSON-file containing the settings stipulated above.

Following files are modified or used as a reference to create files by CoCoNuT, and must be included in the original
OpenFOAM-directory:

- *`system/controlDict`* with compulsory arguments:

|                 key | value                                                                               |
|--------------------:|-------------------------------------------------------------------------------------|
|      `writeControl` | `timeStep`                                                                          |
|     `writeInterval` | required write interval                                                             |
|       `writeFormat` | `ascii`                                                                             |
|        `timeFormat` | `fixed`                                                                             |
|     `timePrecision` | required time precision based on `delta_t` (used in the names of time step folders) |
| `runTimeModifiable` | `false`                                                                             |
|    `adjustTimeStep` | `no`                                                                                |

- *`constant/dynamicMeshDict`* which contains the settings for OpenFOAM's dynamic motion solver
- *`system/decomposeParDict`* with the necessary decomposition of the fluid domain (if `cores`>1)
- *`0/pointDisplacement`* with all the boundary conditions, including `fixedValue` boundary condition for the FSI
  boundaries. This is used as a template for the *`pointDisplacementTmp`* to supply displacement boundary condition (
  from structural solver) for the FSI-interface.

### Comments

- It is probably best to derive a new case from the directory containing an FSI simulation with OpenFOAM in *`coconut/examples/`* in order to copy its structure.
- For OpenFOAM 8, the applications like `pimpleFoam` need to be adapted to accommodate the communication with the solver wrapper during the FSI-simulation.
  These adapted versions have the same name as the original OpenFOAM-solver but with the prefix `coconut_`.
  If you do not use an OpenFOAM-application which is already converted for operation in CoCoNuT, you will have to convert the application yourself, see also [this brief set of instructions](#adapting-a-new-application-to-be-used-in-coconut--openfoam-8-).
- For OpenFOAM 11, only `foamRun` has been adapted but this solver can work with many solver modules (see [OpenFOAM 11 documentation](https://doc.cfd.direct/openfoam/user-guide-v11/solvers-modules)). 

## Treatment of kinematic pressure and traction
If the OpenFOAM application is intended for incompressible flows, e.g. `pimpleFoam` (v8) or `incompressibleFluid` (v11), it solves the incompressible Navier-Stokes equations.
As a consequence the pressure and traction (wallShearStress) values are kinematic and have the unit m²/s².
In order to couple these solvers to a structural solver, the pressure and traction are required to be expressed in Pa (kg/ms²).

For these solvers, the `density` parameter is required in the JSON file and will be used to calculate the actual values from the kinematic ones.
In case the OpenFOAM application is compressible, this correction is not required and the actual values in Pa are obtained directly.

There is a third possibility. In some cases, the application solves the compressible equations obtaining an actual pressure,
but the used momentumTransportModel is incompressible.
It is the type of momentumTransportModel that determines whether the wallShearStress functionObject returns the kinematic or actual traction.
In these cases, the density is not fixed and cannot simply be multiplied with the obtained kinematic values.
Therefore, a new functionObject, rhoWallShearStress, is available, which is a modified version of wallShearStress
and returns the actual traction, even if the momentumTransportModel is incompressible.
Note that is this new functionObject is version specific, but common for all applications of one version.
Upon compilation of the coconut application this functionObject is included in the binary.
This is the case for, for example, coconut_interFoam and coconut_cavitatingFoam.

The behaviour of CoCoNuT for these different applications is implemented in the solver wrapper itself, namely in the `kinematic_conversion_dict`, located in the files *`vX.py`*, with `X`
the identifier of the OpenFOAM-version (e.g. *`v8.py`* for OpenFOAM 8).
This means that when a new application is added, the behaviour for this new application should be specified in that location.

## Solver coupling convergence

The convergence criterion [solver coupling convergence](../../convergence_criteria/convergence_criteria.md#solver-coupling-convergence) has been implemented for the OpenFOAM solver wrapper.
Instead of comparing the different residuals (U, p, ...) with their respective tolerances before the first iteration, they are checked after the first iteration due to implementation restrictions.
Note that OpenFOAM will always execute a final solver iteration after convergence.
Depending on the settings, this final iteration may have different or no relaxation factors.
Calling the solver twice with the same interface displacement, will not necessarily mean that the second call immediately converges.
This is due to a (weak) nonlinearity in the calculation of the mesh deformation as a result of the explicit non-orthogonality correction in the normal gradient calculation used in the Laplacian equation.
Attempts to use the option moveMeshOuterCorrectors (recalculating mesh motion in every solver iteration) showed that there exist problems with this feature: it disturbs the convergence checks and, when used in parallel, even causes the code to crash.

## Version specific documentation

### v8 (OpenFOAM 8)

Base version.
In this OpenFOAM version some changes were made with respect to older versions, for example, the name of some commands and classes are different.
Moreover, there are some changes in the case setup. The following list provides some but is not all extensive:

- The interpolation point for `arc` definitions using `blockMesh` has to be on the arc and not just on its 'elongation'
- Using a macro expansion, such as `$R`, in combination with `#calc`, where an operator follows the macro expansion,
  requires using protection: use `#calc "${R}/2` instead of `#calc "$R/2`
- The file _`constant/turbulenceProperties`_ has been renamed to _`constant/momentumTransport`_
- The keyword `residualControl` in the `PIMPLE` dictionary is now called `outerCorrectorResidualControl`. The
  keyword `residualControl` still exists, but has a different meaning and is used as in `SIMPLE`.
- The file _`system/fvSolution`_ requires dictionaries `pcorr` and `pcorrFinal`, if the keyword `correctPhi` is `true`
  in _`system/fvSolution`_, which it is by default. The tolerance settings should be sufficiently low (similar to `p`).

### v11 (OpenFOAM 11)

In this OpenFOAM version some changes were made with respect to OpenFOAM 8.
The most prominent is the replacement of separate applications like `pimpleFoam` by a single solver `foamRun` which works with solver modules.
The required changes for the OpenFOAM wrapper are visible in the version specific solver wrapper Python file *`v11.py`*.

Moreover, there are some changes in the case setup. The following list provides some but is not all extensive:

- The application in _`system/controlDict`_ is now foamRun and an additional keyword solver is required, such as `incompressibleFluid`.
- The keyword `version` can be omitted from the header dictionary `FoamFile`.
- The file _`constant/transportProperties`_ has been renamed to _`constant/physicalProperties`_.
- The dynamicMeshDict has been restructured, defining a `fvMeshMover` and removing the `dynamicFvMesh` class, see also the [OpenFOAM documentation](https://cfd.direct/openfoam/free-software/dynamic-meshes/).

## Details on the OpenFOAM solver wrapper

### Communication with OpenFOAM

As with most other solver wrappers, the communication between the OpenFOAM code and the solver wrapper occurs 
through the use of (empty) files with extension `.coco`.
The OpenFOAM code checks for these message files in an infinite loop.
When such file is detected, for example _`continue.coco`_, 
the OpenFOAM code performs a certain action and creates a file _`continue_ready.coco`_.
The solver wrapper is paused until it detects such a corresponding file.

### Clean-up after unexpected stop

The aforementioned messaging procedure implies that OpenFOAM is constantly running during execution of the CoCoNuT-simulation.
Exiting of the loop and closing of the program only occurs when the _`stop.coco`_ is received.
If an unexpected crash occurs, it can occasionally occur that the OpenFOAM processes are not ended properly 
and that the user should take care to kill that OpenFOAM-loop manually (using `kill` or `pkill` in the Linux-terminal, e.g. `pkill coconut_*`).

### The `solve_solution_step` method

This method is the core of the simulation passing on the interface displacement to OpenFOAM and receiving the resulting loads (pressure and traction).
The interface displacement is converted into an OpenFOAM-readable format (with the method `write_node_input`),
by storing it in a `pointDisplacementTmp` field, which is read in by OpenFOAM in every iteration 
(this required some adaptation of the solver, see [next section](#adapting-a-new-application-to-be-used-in-coconut--openfoam-8-)). 
Subsequently, the mesh is updated in OpenFOAM and the flow equations are solved.
The dynamic mesh motion is handled by OpenFOAM itself.
Finally, the method `read_node_output` is called in the solver wrapper, which reads the interface loads from the directory *`postProcessing`* (more precisely from the subdirectories *`coconut_<boundary_name>`*).

## Adapting a new application to be used in CoCoNuT (OpenFOAM 8)

For OpenFOAM 8, the applications like `pimpleFoam` need to be adapted to accommodate the communications with the solver wrapper during the FSI-simulation.
These adapted versions have the same name as the original OpenFOAM-solver but with the prefix `coconut_`.
If you do not use an OpenFOAM-application which is already converted for operation in CoCoNuT, you will have to convert the application yourself, see also []().
This can be done in a rather straightforward way by taking a look at already implemented application, for example `coconut_pimpleFoam`.
In brief, the following steps should be undertaken:

- Some additional `include`-statements are needed: _`fsiDisplacement.H`_, _`waitForSync.H`_ and _`<unistd.h>`_.
  Also include _`readCoconutControls.H`_ to create control variables used by CoCoNuT.
- Except for these initial `include`-statements, the entire solver code should be put in an *infinite* loop that starts
  with `while (true)`. The code will stay in this loop and check with several conditional statements whether the solver wrapper in CoCoNuT has
  sent a message to the OpenFOAM-solver.
  These messages are sent by creating an empty file with a specific name in the   OpenFOAM-directory.
  The following file names should be checked by the OpenFOAM-solver: *`next.coco`*, *`continue.coco`*, *`save.coco`*, *`stop.coco`*.
  To allow a pause of 1 ms between each loop the command usleep(1000) is used (which requires the line: `#include <unistd.h>` before the loop).
- Once such a message is received, the OpenFOAM code first executes a command to sync all processors (only important in a parallel run), for example `waitForSync("next")`.
  Details on this function are given [here](#the-waitforsync-command).
- If the file *`next.coco`* exists, a new time step is started.
  Due to the fact that the infinite loop is driven by `while(true)` and not the default OpenFOAM-expression `while(runTime.run())`, 
  the statement `runTime.run();` has to be added. 
  This creates, initializes and calls the functionObjects in the *`controlDict`* file.
  Furthermore, as in the original application, the runTime-object should be incremented, increasing the time step.
- If the file *`continue.coco`* exists, the interface has to be moved first, using the following lines:
  ``` cpp
  forAll(boundaryNames, s)
  {
      word boundaryName = boundaryNames[s];
      ApplyFSIPointDisplacement(mesh, boundaryName);
  }
  ```
  Thereafter, the flow equations need to be solved. This `if`-statement consequently contains
  most of the original solver definition, in which the flow equations are called in the same order as in the original
  CFD solver. Additionally, at the end of the PIMPLE corrector loop, the coupling convergence is checked, see [solver coupling convergence](#solver-coupling-convergence).
  Finally, at the end of this block the loads are written by calling the CoCoNuT functionObjects by including _`executeCoconutFunctionObjects.H`_.
- If the file *`save.coco`* exists, it is checked whether the flow fields should be stored in corresponding files according to the user-defined save interval by calling _`runTime.write();`_.
- If the file *`stop.coco`* exists, a `break`-statement should end the infinite loop and the OpenFOAM program terminates.

The solver wrapper will automatically compile the adapted applications.
This can also be done the default OpenFOAM-compilation method (loading the OpenFOAM-module and using `wmake` in the directory *`solver_wrapper/coconut_<solver_name>`*).
To be able to compile the new application, the following files need to adapted.

- Update the file _`Make/files`_ to have the correct application name and change the location of the executable to `$(FOAM_USER_APPBIN)`:
  ``` cpp
  coconut_pimpleFoam.C

  EXE = $(FOAM_USER_APPBIN)/coconut_pimpleFoam
  ```
- Update the file _`Make/options`_ such that the CoCoNuT header files are found by adding the parent directory and its parent directory to the included files:
  ``` cpp
  EXE_INC = \
    -I.. \
    -I../.. \
    -I<INSERT OTHER LOCATIONS>
  ```

### The waitForSync command

This command, for example `waitForSync("next")`, is called at the start of every message block and is required to avoid that one processor goes through the block and already removes the message file before all others have seen it.
This is achieved by gathering a label on all processors, effectively syncing them.
After they are synced, the message file is removed and the corresponding _ready_ message is written, for example `next_ready.coco`.
Thereafter, they are synced once more.
This avoids that while removing the file, one processor has left and again entered the same block.
Note that the solver wrapper already receives the return message, for example `next_ready.coco`, before the completion of the block.
This results in a speed-up, but requires the solver wrapper to check if, for example, the loads have been fully written before reading them.

## Disclaimer

This offering is not approved or endorsed by OpenCFD Limited, producer and distributor of the OpenFOAM software via www.openfoam.com, and owner of the OPENFOAM® and OpenCFD® trademarks,
nor by OpenFOAM Foundation Limited, producer and distributor of the OpenFOAM software via www.openfoam.org.