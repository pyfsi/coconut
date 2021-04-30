# Turek benchmark with Fluent2D and Abaqus2D

This example calculates the well-known Turek benchmark, which consist laminar incompressible flow around a flexible flag attached to a rigid cylinder in a 2D channel.
The material parameters and boundary conditions correspond to the FSI2 setup, as detailed by Turek and Hron [[1](#1)].
The challenging aspect of this benchmark lies in the large self-indcuced oscillation of the flag.
The used solvers are Fluent and Abaqus.

Due to the time required for the oscillations to enter a periodic regime, this test case takes a long time to run.

The figure below shows the result.

In the FSI2 setup a parabolic velocity profile is subscribed at the inlet with an average velocity of 1 m/s, ramped up slowly in time.  
The fluid paramters are: 

-   density: 1000 kg/m³
-   dynamic viscosity: 1 Pa$\cdot$s

The flag is consist of a linear elastic material with the follwing properties:

-   density: 10$^4$ kg/m³
-   modulus of elasticity: 1.4 10$^6$ Pa
-   poisson's ratio: 0.4

The total simulated time is 35 s in time steps of 0.002 s.

## Coupling algorithm

The coupling technique used is the *interface quasi-Newton algorithm with an approximation for the inverse of the Jacobian from a least-squares model* (IQNI-LS).
The reuse parameter `q` is set to 10, which means that data from the last ten time steps will be used to stabilize and accelerate the convergence.

## Predictor

The initial guess in every time step is done using the quadratic predictor.

## Convergence criterion

Two convergence criteria have been specified:

-   The number of iterations in every time step is larger than 15.
-   The residual norm on the displacement is a factor $10^{-3}$ lower than the initial value.
 
When either criterion is satisfied the simulation stops.

## Solvers

Fluent is used as flow solver.
The mesh is provided. It is triangular and rather coarse to reduce the calculation time.
A script to regenerate it using Gambit is included. This script allows to change the resolution and geometrical parameters.

The structure solver is Abaqus.
The Abaqus case is built when setting up the case starting from the file *`mesh_turek.inp`* containing nodes and elements. 
This is done by running Abaqus with the *`makeInp.py`* Python script to set all parameters, such as surface definitions, material parameters, boundary conditions and time step information.
The result of the setup is a completed input file *`case_turek.inp`*.
The Abaqus element type used is CPE8R.

To exchange information between the solvers on the fluid-structure interface, the use of mappers is required.
In the structure solver wrapper, a linear interpolation mapper is used to interpolate in the x- and y-direction from and to the coupled solver.
When the case is run, a bounding box warning for the right most face of the flag appears.
This is normal.
The reason is that the distance between the load points in Abaqus and the face centers in Fluent is large compared to the overall dimension of this face, 
due to the coarse resolution: only 3 cells in the flow mesh and 3 elements in structure mesh.

## References
<a id="1">[1]</a> 
[Turek S., Hron J., "Proposal for Numerical Benchmarking of Fluid-Structure Interaction between an Elastic Object and Laminar Incompressible Flow", Bungartz HJ., Schäfer M. (eds) Fluid-Structure Interaction. Lecture Notes in Computational Science and Engineering, vol. 53. Springer, Berlin, Heidelberg, 2006.](https://doi.org/10.1007/3-540-34596-5_15)