# Turek benchmark with Fluent2D and Abaqus2D - Steady

This example calculates the well-known Turek benchmark, which consist of laminar incompressible flow around a flexible flag attached to a rigid cylinder in a 2D channel.
The material parameters and boundary conditions correspond to the FSI1 setup, as detailed by Turek and Hron [[1](#1)].
This is a steady simulation.
The used solvers are Fluent and Abaqus.
A script _`evaluate_benchmark.py`_ is provided to compare the results with the benchmark results in [[1](#1)].

The figure below shows the resulting velocity contour plot (with Fluent).
![FSI1](images/turek_fsi1_velocity.png "Velocity contour plot of FSI1 setup produced with Fluent")

In the FSI1 set up a parabolic velocity profile is subscribed at the inlet with an average velocity of 0.2 m/s.
The flow is initialized with a velocity in the x-direction equal to 0.2 m/s.  
The fluid parameters are:

-   density: 1000 kg/m³
-   dynamic viscosity: 1 Pa$\cdot$s

The flag is consist of a linear elastic material with the following properties:

-   density: 1000 kg/m³
-   modulus of elasticity: 1.4 10$^6$ Pa
-   poisson's ratio: 0.4

## General setting
Although the calculation is steady, a `delta_t` is still required. Its value is arbitrary and 0 is used.
`timestep_start` is required as well and is set to 0. The `number_of_timesteps` is set to 1.

## Coupling algorithm

The coupling technique used is the *interface quasi-Newton algorithm with an approximation for the inverse of the Jacobian from a least-squares model* (IQN-ILS).
Note that the reuse parameter `q` is not used as only one time step is calculated.

## Predictor

A predictor is still required, but not used as only one time step is calculated.

## Convergence criterion

Two convergence criteria have been specified:

-   The number of iterations in every time step is larger than 15.
-   The residual norm of the displacement is a factor $10^{-3}$ lower than the initial value.
 
When either criterion is satisfied the simulation stops.

## Solvers

Fluent is used as flow solver.
The provided mesh is triangular.
A script to regenerate it using Gambit is included. This script allows to change the resolution and geometrical parameters.

The structural solver is Abaqus.
The Abaqus case is built when setting up the case starting from the file *`mesh_turek.inp`* containing nodes and elements. 
This is done by running Abaqus with the *`make_inp.py`* Python script to set all parameters, such as surface definitions, material parameters, boundary conditions and time step information.
The result of the setup is a completed input file *`case_turek.inp`*.
The Abaqus element type used is CPE8R.

To exchange information between the solvers on the fluid-structure interface, the use of mappers is required.
In the structural solver wrapper, a linear interpolation mapper is used to interpolate in the x- and y-direction from and to the coupled solver.

## References
<a id="1">[1]</a> 
[Turek S., Hron J., "Proposal for Numerical Benchmarking of Fluid-Structure Interaction between an Elastic Object and Laminar Incompressible Flow", Bungartz HJ., Schäfer M. (eds) Fluid-Structure Interaction. Lecture Notes in Computational Science and Engineering, vol. 53. Springer, Berlin, Heidelberg, 2006.](https://doi.org/10.1007/3-540-34596-5_15)