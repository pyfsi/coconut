# Tube case with TubeFlow and TubeStructure

This example calculates the flow inside and the deformation and stresses of a straight flexible tube, where a sinusoidal velocity is specified at the inlet.
This done by using Python solvers TubeFlow and TubeRingmodel.

## Coupling algorithm

The coupling technique used is the *interface quasi-Newton algorithm with an approximation for the inverse of the Jacobian from a least-squares model* (IQNI-LS).

## Predictor

The initial guess in every time step is done using the linear predictor.

## Convergence criterion

Two convergence criteria have been specified:

-   The number of iterations in every time step is larger than 15.
-   The residual norm on the displacement is a factor $10^{-6}$ lower than the initial value.
 
When either criterion is satisfied the simulation stops.

## Solvers

The flow solver is the Python solver TubeFlow, which implements a 1D model of the flow inside the tube,
with 100 cells on the fluid-structure interface. 
The parameters for this model are specified in the setup folder by the file `solver_parameter.json`.
The (radial) displacements are applied on the cell centers.
The loads, in fact are only pressure for this 1D case, are calculated in the cell centers as well.
The axial direction is along the z-axis,
the radial direction along the y-axis.

The structure solver is the Python solver TubeRingmodel, which implements a 1D model of the tube wall,
with 100 elements on the fluid-structure interface.
It differs from the Python solver TubeStructure as the no inertia is considered.
The tube is regarded as consisting out of 100 independent rings.
The parameters for this model are specified in the setup folder by the file `solver_parameter.json`.
The loads, which in fact only pressure for this 1D case, are applied on the cell centers.
The displacements are calculated in the cell centers as well.
Only the radial displacement is different from zero.
The axial direction is along the z-axis,
the radial direction along the y-axis.

The coordinate frames of both solvers are the same so there is no need for a permutation mapper.
However, the discretization of both solvers differ.
To account for the difference of the points where loads and displacements are applied or calculated, the use of interpolation mappers is required.
Therefore, a linear interpolation mapper is introduced in the structure solver to interpolate in the x-direction.

Additionally a parameter file `parameters_conformal.json` is also provided.
This parameter file performs no interpolation.
It should be verified that the number of cells `m` is the same in both solvers before using it and is merely provided as a theoretical example as it will have close to no practical use.