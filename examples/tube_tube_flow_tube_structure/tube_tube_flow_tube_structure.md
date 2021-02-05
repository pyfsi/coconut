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
The parameters for this model are specified in the setup folder by the file `solver_parameter.json`.
The (radial) displacements are applied on the cell centers.
The loads, in fact are only pressure for this 1D case, are calculated in the cell centers as well.
The axial direction is along the z-axis,
the radial direction along the y-axis.

The structure solver is the Python solver TubeStructure, which implements a 1D model of the tube wall,
with 100 elements on the fluid-structure interface.
The parameters for this model are specified in the setup folder by the file `solver_parameter.json`.
The loads, which in fact only pressure for this 1D case, are applied on the cell centers.
The displacements are calculated in the cell centers as well.
Only the radial displacement is different from zero.
The axial direction is along the z-axis,
the radial direction along the y-axis.

The coordinate frames of both solvers are the same so there is no need for a permutation mapper.
As both solver has the 100 cells on the fluid-structure interface, no interpolation is required
and the parameter file `project_parameters_conformal.json` can be used.

A parameter file `project_parameters_mapped.json` is also provided,
which can be used i the number of cells `m` is set differently in the two solvers.
Then, a linear interpolation mapper is introduced in the structure solver to interpolate in the x-direction.