# Tube case with Fluent3D and KratosStructure3D

This example calculates the flow inside and the deformation and stresses of a straight flexible tube, where a pressure pulse is applied at the inlet.
This done by using Fluent and StructuralMechanicsApplication of Kratos, both with a fully 3D case.

## Coupling algorithm

The coupling technique used is the *interface quasi-Newton algorithm with an approximation for the inverse of the Jacobian from a least-squares model* (IQNI-LS).

## Predictor

The initial guess in every time step is done using the linear predictor.

## Convergence criterion

Two convergence criteria have been specified:
 - The number of iterations in every time step is larger than 20.
 - The residual norm on the displacement is a factor $10^{-3}$ lower than the initial value.
 
When either criterion is satisfied the simulation stops.
 
## Solvers

The flow solver is Fluent, used to solve a fully 3D tube,
with 48 cells on the fluid-structure interface in the length-wise direction and 8 in the circumferential direction.
When setting up the case, the mesh is build based on the file `mesh.jou` using Gambit.
The displacements are applied in the nodes. 
In contrast, the loads (pressure and traction) are calculated in the cell centers.
The axial direction is along the x-axis,
the radial direction along the y-axis.

The structure solver is Kratos, used to solve a fully 3D tube,
with 24 elements on the fluid-structure interface in the length-wise direction and 8 in the circumferential direction. 
The Kratos element type used is ShellThickElementCorotational3D4N. These are 4-node shell elements.
The displacement and loads are calculated/applied on the nodes. There are 200 nodes on the fluid-structure interface.
The axial direction is along the x-axis,

The coordinate frames of both solvers are the same so there is no need for a permutation mapper.
On the other hand, difference in the discretization between the fluid and the structure mesh requires the use of interpolation mappers.
Therefore, a radial basis mapper is used to interpolate in the x-, y- and z-direction.