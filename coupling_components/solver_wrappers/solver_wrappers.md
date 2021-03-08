# Solver-wrappers

The goal of solver-wrappers is to provide communication with the solvers. This  means that a solver-wrapper must implement a way to communicate input to the solver, run the solver with the provided input and obtain the solution from the solver. 

As each solver is different, the solver-wrappers are highly customized too. Nevertheless, they all inherit from the `Component` class and must all implement the following methods: `initialize`, `initialize_solution_step`, `solve_solution_step`, `finalize_solution_step`, `finalize`, `get_interface_input` and `get_interface_output`. 

The `solve_solution_step` method is called each coupling iteration with input data on the interface, which it provides to the solver. After the solver has obtained a solution, it reads this solution and returns it as output data on the interface.

## Implementations

There are currently 2 solver-wrappers for flow solvers implemented:

-   [ANSYS Fluent](fluent/fluent.md) (2019R1, 2019R2, 2019R3, 2020R1)
-   [OpenFOAM](openfoam/openfoam.md) (4.1)

There are currently 2 solver-wrappers for structural solvers implemented:

-   [Abaqus](abaqus/abaqus.md) (6.14)
-   [Kratos Multiphysics](kratos/kratos.md) (6.0)
