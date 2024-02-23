# Breaking dam with Fluent2D and KratosStructure2D

This example simulates a dam-break on a top-fixed flexible gate.
The liquid, initially located in a free-surface tank, imposes a variable pressure on the bending flexible gate.
This case has been examined numerically and experimentally by Antoci et al. [[1](#1)].
Here, this 2D FSI calculation is performed with Fluent and Kratos Multiphysics.

For details about geometry, parameters and boundary conditions refer to the example with [Fluent and Abaqus](../fluent2d_abaqus2d/breaking_dam_fluent2d_abaqus2d.md).

## Coupling algorithm

The coupling technique used is IQN-ILSM with reuse of $q = 50$ time steps.
These are used to stabilize and accelerate the convergence.

## Predictor

The initial guess at each time step is made using the linear predictor.

## Convergence criterion

Two convergence criteria have been specified:

-   The number of iterations in every time step is larger than 20.
-   The residual norm of the displacement is a factor $10^{3}$ lower than the initial value.

The simulation stops as soon as one of the two criteria is met.

## Solvers

Fluent is used as flow solver.
The provided mesh is triangular. When the gate bends, remeshing is performed to preserve its quality.
A script to regenerate it using Gambit is included. This script allows the resolution and geometric parameters to be
changed.

The structural solver is Kratos Multiphysics. First-order, quadrilateral, plane strain elements are used.

To exchange information between the solvers on the fluid-structure interface, the use of mappers is required.
In the structural solver wrapper, a linear interpolation mapper is used to interpolate in the x- and y-direction from and to the coupled solver.

## References
<a id="1">[1]</a>
[Antoci C., Gallatie M. and Sibilla S., "Numerical simulation of fluid-structure interaction by SPH", Computers & Structures, vol. 85, no. 11, pp. 879-890, 2007.](https://doi.org/10.1016/j.compstruc.2007.01.002)