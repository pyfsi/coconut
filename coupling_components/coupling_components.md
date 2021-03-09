# Coupling components

The coupling components are the basic building blocks of the CoCoNuT coupling tool. They are defined in the subdirectories of the *`$COCO/coconut/coupling_components`* directory.
There are five types: 

-   [*Convergence criteria*](convergence_criteria/convergence_criteria.md) which determine when the calculation within a time step has converged. Subdirectory: *`convergence_criteria`*.
-   [*Coupled solvers*](coupled_solvers/coupled_solvers.md) (completed by *models*, subdirectory *`coupled_solvers/models`*) which perform the actual coupling, implementing a coupling algorithm. Subdirectory *`coupled/solvers`*.
-   [*Mappers*](mappers/mappers.md) which map from one interface discretization to another, i.e. interpolation between non-conformal meshes. Subdirectory *`mappers`*.
-   [*Predictors*](predictors/predictors.md) which provide an initial guess at the start of a new time step. Subdirectory *`predictors`*.
-   [*Solver-wrappers*](solver_wrappers/solver_wrappers.md) which provide communication which the actual solvers. Subdirectory *`solver_wrappers`*.

The idea behind these components is modularity. For example, changing a solver wrapper or creating a new one can be done without having to adapt any other components.
This allows for high degree of flexibility, a short learning curve and a limited development effort.
Moreover, components can be used multiple times without the need for copying code.
Detailed information on these components can be found in the specific documentation.

All these coupling components inherit from a same superclass called `Component`, defined in *`$COCO/coconut/coupling_components/component.py`*. 
In this class some methods are implemented, which control the flow within CoCoNuT. 
For every coupling component, there are `initialize` and `finalize`, which are called at the start and end of a calculation
and there are `initialize_solution_step` and `finalize_solution_step`, which are called at the start and end of each time step.
If needed these methods are overwritten to perform a specific action.
For example in the `finalize` method of a `solver_wrapper`, the termination of the solver software could be implemented.

A schematic of the relation between the coupling components for a basic calculation is given in the following figure.
![](images/coupling_components.png "Schematic of relation between coupling components")

These coupling components have to communicate with each other.
This is done through the use of *interfaces*, indicated with arrows on the figure. 
For these `Interface` objects (implemented in *`$COCO/coconut/data_structure/interfaces.py`*) containing the (discretized) solution data on the FSI-interface and references to among others, the coordinates of the discretized interface. The implementation of this `Interface` class is explained in more detail in [the documentation about the data structure](../data_structure/data_structure.md).

## Start of the calculation

The main coupling component in which all other coupling components are instantiated is the coupled solver.
The coupled solver itself is created in the `Analysis` class (*`$COCO/coconut/analysis.py`*, which is the starting point of the CoCoNuT calculation.
Upon the start of CoCoNuT, an instance of `Analysis` is made and its method `run` is executed.
The coupled solver keeps track of all `Components` and runs the methods `initialize`, `finalize`, `initialize_solution_step` and `finalize_solution_step`,
when its respective methods are executed.

## Tools

Some code to perform specific tasks, like printing with a certain layout or performing a time measurement is useful throughout the code.
These functionalities are grouped in the file *`$COCO/coconut/tools.py`*.
It suffices to import the file to make uses of its functions.
