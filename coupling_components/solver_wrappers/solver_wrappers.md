# Solver wrappers


The goal of solver wrappers is to provide communication with the solvers. This  means that a solver wrapper must implement a way to communicate input to the solver (at the fluid-structure interface), run the solver with the provided input and obtain the solution from the solver (at the fluid-structure interface).

**Important**: To avoid conflicts, each solver is run in its own environment. These environments are established on runtime, but the actions required to do so can differ between systems. Therefore it is necessary to first check (and if needed adapt) the file *`coconut/solver_modules.py`* after installing or updating CoCoNuT, as described in [the front documentation page](../../README.md#checking-the-solver-modules).

As each solver is different, the solver wrappers are highly customized too. Nevertheless, they all inherit from the `Component` class and must all implement the following methods: `initialize`, `initialize_solution_step`, `solve_solution_step`, `finalize_solution_step`, `finalize`, `get_interface_input` and `get_interface_output`.

Upon instantiation of the solver wrapper object, the solver wrapper has to create a `Model` containing one or more `ModelParts` which correspond to (a part of) the fluid-structure interface. The interface coordinates at time step 0 should be stored in the `ModelPart`, even for a simulation that is restarted from a later time step, in order to have consistent [mapping](../mappers/mappers.md) over restarts.

The `solve_solution_step` method is called in each coupling iteration by the coupled solver, which provides input data stored in an `Interface` object. The solver wrapper extracts the data and supplies it to the actual solver, starts the calculation and reads the output data when the solver has finished. The solver wrapper then returns this data to the coupled solver as an `Interface` object.


## Available solver wrappers

There are currently two solver wrappers for computational fluid dynamics (CFD) packages:

-   [ANSYS Fluent](fluent/fluent.md) (2019R1, 2019R2, 2019R3, 2020R1)
-   [OpenFOAM](openfoam/openfoam.md) (4.1)

There are currently also two solver wrappers for computational structural mechanics (CSM) packages :

-   [Abaqus](abaqus/abaqus.md) (6.14)
-   [Kratos Multiphysics](kratos_structure/kratos_structure.md) (6.0)

CoCoNuT also implements several 1D Python-based solver wrappers to provide a fast and easy way to test new algorithms and implementations. These are:

-   [Tube flow solver](python/python.md#tube-flow-solver)
-   [Tube structure solver](python/python.md#tube-structure-solver)
-   [Tube ringmodel solver](python/python.md#tube-ringmodel-solver)