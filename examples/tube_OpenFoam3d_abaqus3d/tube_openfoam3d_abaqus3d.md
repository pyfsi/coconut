# Tube case with OpenFOAM3D and Abaqus3D
!! Should still be adapted to OpenFOAM, was copied from Fluent!!
This example calculates the flow inside and the deformation and stresses of a straight flexible tube, where a pressure pulse is applied at the inlet.
This done by using OpenFOAM and Abaqus, both with a fully 3D case.

## Coupling algorithm

The coupling technique used is the *interface quasi-Newton algorithm with an approximation for the inverse of the Jacobian from a least-squares model* (IQNI-LS).

## Predictors

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

The structure solver is Abaqus, used to solve a fully 3D tube,
with 12 elements on the fluid-structure interface in the length-wise direction and 8 in the circumferential direction.
The Abaqus case is not build when setting up the case, but is provided as the file `Base.inp`. 
The Abaqus element type used is S8R. These are 8-node shell elements. The R refers to reduced integration.
See Abaqus documentation for more information. 
The loads are applied on the faces in 9 points per element, which means on 864 load points in total. 
The displacement is calculated in the nodes. There are 304 nodes on the fluid-structure interface.
The axial direction is along the x-axis,
the radial direction along the y-axis.

The coordinate frames of both solvers are the same so there is no need for a permutation mapper.
In contrast, the difference of the points where loads and displacements are applied or calculated,
require the use of interpolation mappers.
Therefore, a radial basis mapper is introduced in the structure solver to interpolate in the x-, y- and z-direction.