# Example case to test individual solvers

This example shows how the separate solvers or their cases in the respective working directories can be tested individually
using `test_single_solver` as coupling component. For more information refer to [coupled_solvers](../coupling_components/coupled_solvers/coupled_solvers.md).

## Coupling algorithm

The `type` is set to `coupled_solvers.test_single_solver` and an additional required dictionary `test_settings` is added.
All other parameters in the JSON file, including the `settings` dictionary, can be set to the values that will be used in the actual coupled simulation.
The `settings` dictionary is used to look up `delta_t`, `timestep_start`, `save_results` and `name` if not provided in `test_settings`.
The other dictionaries are not used: no `predictor`,` convergence_criterion` or `mapper` are used.

## Solvers

If you run this case as is, the Fluent case and solver will be tested, as the `solver_index` is set to 0 in `test_settings`.
In that case the Abaqus settings won't be used.
The Abaqus case can easily be tested by changing this value to 1.
Then the Fluent settings will not be used.
The index refers to the index of the respective solver in the `solver_wrappers` list in the JSON file.

## Dummy solver

Even though the presence `dummy_solver.py` with a test class is not strictly required, it can be very valuable because it allows to have a non-zero input.
The file included in this case provides an example. It contains, among others, `SimpleTest`, `TransientTest` and `InterpolatedData` to illustrate how variables can be defined based on undeformed coordinates (`x`, `y`, `z`) and time step (`n`).
Note that the number of classes defined is not restricted.
Also note that the `__init__()` mehtod can be used to avoid repeating the same calculations multiple times.
The name of the class that one wants to use should be specified in the JSON file. 
In this example the `TransientTest` is used.
Each class contains function definitions with a pre-formatted name `calculate_VARIABLE(x,y,z,n)`, with `<variable>` being the variable(s) required by the tested solver, e.g. `DISPLACEMENT`, `PRESSURE` or `TRACTION`.
How these variables are defined inside these methods, is up to the user.
However, the methods need to return the right format: a 3-element list of floats for vector variables and a single float for scalar variables.