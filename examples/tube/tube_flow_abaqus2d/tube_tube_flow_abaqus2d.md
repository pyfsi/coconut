# Tube case with TubeFlow and Abaqus2D

This example calculates the flow inside and the deformation and stresses of a straight flexible tube, where a pressure pulse is applied at the inlet.
This done by using Python solver TubeFlow and Abaqus (axisymmetric).

## Coupling algorithm

The coupling technique used is the *interface quasi-Newton algorithm with an approximation for the inverse of the Jacobian from a least-squares model* (IQN-ILS).

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
The axial direction is along the z-axis,
the radial direction along the y-axis.

The structural solver is Abaqus, used to solve an axisymmetric representation of the tube,
with 50 elements on the fluid-structure interface.
The Abaqus case is built when setting up the case starting from the file *`mesh_tube2d.inp`* containing nodes and elements. 
This is done by running Abaqus with the *`make_inp.py`* Python script to set all parameters, such as surface definitions, material parameters, boundary conditions and time step information.
The result of the setup is a completed input file *`case_tube2d.inp`*.The Abaqus element type used is CAX8RH. These are continuum elements for axisymmetric calculations, for stress and displacement without twist. 
They are: 8-node biquadratic, reduced integration, hybrid with linear pressure. 
See the [Abaqus documentation](http://130.149.89.49:2080/v6.14/books/usb/default.htm?startat=book01.html#usb) for more information. 
The loads are applied on the faces in three points per element, which means on 150 load points in total. 
The displacement is calculated in the nodes. There are 101 nodes on the fluid-structure interface.
The axial direction is along the y-axis, the radial direction along the x-axis.

The difference in reference frames and number of cells on the fluid-structure interface requires the use of mappers.
In the structural solver wrapper, a permutation mapper is introduced to match the coordinate frames.
The new x-axis corresponds to y-axis of the incoming loads the new y-axis corresponds to the z-axis and the new z-axis to the x-axis.
Thereafter, a linear interpolation mapper is used to interpolate in the x- and y-direction from the 100 cell centers of TubeFlow to the 150 load points in Abaqus.
For the output the same is done in the opposite order: first interpolating and then swapping the axes.
So the new x-axis is the z-axis of the displacements calculated by Abaqus, the new y-axis corresponds to the x-axis and the new z-axis to the y-axis.