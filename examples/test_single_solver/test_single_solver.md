# Example case to test individual solvers

This example shows how the separate solvers or their cases in the respecive working directories can be tested individually
using `test_single_solver` as coupling component (documentation on this can be found under `coupling_components/coupled_sovers`.

## Coupling algorithm

The `type` is set to `test_single_solver` and an additional required dictionary `test_settings` is added. All other parameters
in the JSON file can be set to the values that will be used in the actual coupled simulation.

## Predictor

The initial guess in every time step is done using the linear predictor. 

## Convergence criterion

The convergence criteria are set as they would be used in the fully coupled simulation. However, as a dummy solver 
imposes a certain motion or certain forces on the `interface_input` of the tested solver wrapper regardless of the output
of the solver wrapper, these settings do not matter in the test environment. 

## Solvers

If you run this case as it is, the Fluent case and solver will be tested, as the `solver_index` is set to 0 in `test_settings`.
The abaqus case can easily be tested by changing this value to 1.

## Dummy solver

Even though the presence of `dummy_solver.py` is not strictly required, it can be valuable to have a nonzero input. The
file included in this case provides an example.