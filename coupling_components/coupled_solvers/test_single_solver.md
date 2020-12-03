# Test single solver

This part of the documentation describes how the `test_single_solver` works and can be used to test new cases and settings.
The idea behind this component is that only one of the two solvers is tested, while the other one is replaced by a dummy.

## Settings

The JSON settings for this `test_single_solver` are largely the same as in [`coupled_solvers`](coupled_solvers.md), namely
the dictionaries `type`, `settings`, `predictor`, `convergence_criterion` and `solver_wrappers`. The entry for `type` is
obviously `test_single_solver`, the `settings` that one wants to use in the actual simulation, can be kept. Additionally 
to these directories, a mandatory directory `test_settings` is to be defined as well. An example case can be found 
[here](../../examples/test_single_solver/test_single_solver.md). The possibilities for the `test_settings` directory 
are shown below.

parameter|type|description
---:|:---:|---
`delta_t`|double|(optional) Time step size to be used in the test. Is optional as long as this value is defined in the `settings` dictionary. If a different value is defined in both dictionaries, the one defined in `test_settings` is chosen.
`solver_index`|int|Has a value 0 or 1 and indicates the solver that one wants to test. 0 indicates the first solver wrapper that appears in the JSON-file, 1 the second wrapper.
`test_class`|string|(optional) Refers to the class to use in the `dummy_solver.py` file in your case (for more information, see [examples](../../examples/test_single_solver/test_single_solver.md)).
`timestep_start`|int|(optional) Time step to start from. In this test environment, this is set to 0 by default.

Note that `test_settings` overwrites some parameters defined in `settings`. As such, these dictionaries are completely 
independent. This allows the user to set up a normal case but instead of running a full simulation, a test can easily be 
set up by only adding this `test_settings` dictionary and setting the `type` of `coupled_solvers` to `test_single_solver`.

An attentive reader will notice that in the [example case](../../examples/test_single_solver/test_single_solver.md), 
an additional parameter `check_mapping` is defined. This allows to test the mappers as well but this feature is currently 
not implemented and as such, this value is not used.

In your case directory where the simulation is run, a python file `dummy_solver.py` should be present that contains
at least one `test_class`. These classes allow to define forces or displacement on your structure and acts as "second
solver" to the problem. If this file is not present or `test_class` is not defined in `test_settings`, zero input will
be used. If `dummy_solver.py` does not contain the `test_class` defined in `test_settings`, an error will occur.

## Algorithm

The initialization starts with checking whether or not `dummy_solver.py` is present. If not, a boolean `bool_test` is set
to `False` and a message is printed that the input will be used. Then, the presence of `test_settings` is checked and
the parameters are set. Note once more that the parameters in `test_settings` have priority over these defined in `settings`.
Furthermore, only `delta_t` is of interest. Other `settings` are ignored, `timestep_start` is set to 0 and `save_results` 
to `False`.

As a next step, `delta_t` and `timestep_start` are added to the `settings` of the solver wrapper to be tested. The working
directory is copied to a new directory with the same name and a suffix `_testX` with `X` an integer starting from 0.  
Previous test working directories are not overwritten. The initialization finishes with informing the user that either
the `test_class` is used to provide input or, in absence of this parameter, a zero input will be used.

The `SolveSolutionStep` is fairly straightforward: if the boolean `bool_test` is true and a `test_class` was defined,
the `SolutionStepValue` of the nodes on the `interface_input` of the solver wrapper is set to the values defined in the
`test_class`. The `interface_output` of the solver wrapper is ignored. As such, this test environment can be considered
as one way FSI in some way.
