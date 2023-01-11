# Tube case with TubeFlow and TubeStructure

This example calculates the flow inside and the deformation and stresses of a straight flexible tube, where a pressure pulse is applied at the inlet.
This done by using Python solvers TubeFlow and TubeStructure.

## Coupling algorithm

The coupling technique used is the *interface quasi-Newton algorithm with an approximation for the inverse of the Jacobian from a least-squares model* (IQNI-LS).

## Predictor

The initial guess in every time step is done using the linear predictor.

## Convergence criterion

Two convergence criteria have been specified:

- The number of iterations in every time step is larger than 15.
- The residual norm on the displacement is a factor $10^{-6}$ lower than the initial value.
 
When either criterion is satisfied the simulation stops.

## Solvers

The flow solver is the Python solver TubeFlow, which implements a 1D model of the flow inside the tube,
with 100 cells on the fluid-structure interface. 
The parameters for this model are specified in the setup folder by the file *`solver_parameters_pressure.json`*.
The (radial) displacements are applied on the cell centers.
The loads, in fact only pressure for this 1D case, are calculated in the cell centers as well.
The axial direction is along the z-axis, the radial direction along the y-axis.

The structure solver is the Python solver TubeStructure, which implements a 1D model of the tube wall,
with 100 elements on the fluid-structure interface.
The parameters for this model are specified in the setup folder by the file *`solver_parameters.json`*.
The loads, in fact only pressure for this 1D case, are applied on the cell centers.
The displacements are calculated in the cell centers as well. Only the radial displacement is different from zero.
The axial direction is along the z-axis, the radial direction along the y-axis.

The coordinate frames of both solvers are the same so there is no need for a permutation mapper.
Moreover, as both solvers have 100 cells on the fluid-structure interface, no interpolation is required.
Nonetheless, a linear interpolation mapper for the structure solver to interpolate in the x-direction is included in the parameter file. 
As such, the number of cells `m` can be varied independently in both solvers.