# Turek benchmark with Fluent2D and Abaqus2D

This example calculates the well-known Turek benchmark, which consist of laminar incompressible flow around a flexible flag attached to a rigid cylinder in a 2D channel.
The material parameters and boundary conditions correspond to the FSI2 setup, as detailed by Turek and Hron [[1](#1)].
The challenging aspect of this benchmark lies in the large self-induced oscillation of the flag.
The used solvers are Fluent and Abaqus.
A script _`evaluate_benchmark.py`_ is provided to compare the results with the benchmark results in [[1](#1)].
Furthermore, the script _`animate_beam.py`_ can be used to visualize the interface data throughout the calculation.

Due to the time required for the oscillations to enter a periodic regime, this test case takes a long time to run.

The figure below shows the resulting velocity contour plot (with Fluent).
![FSI2](images/turek_fsi2_velocity.gif "Velocity animation of FSI2 setup produced with Fluent")

In the FSI2 set up a parabolic velocity profile is subscribed at the inlet with an average velocity of 1 m/s, ramped up slowly in time.  
The fluid parameters are: 

-   density: 1000 kg/m³
-   dynamic viscosity: 1 Pa$\cdot$s

The flag is consist of a linear elastic material with the following properties:

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
-   The residual norm of the displacement is a factor $10^{-3}$ lower than the initial value.
 
When either criterion is satisfied the simulation stops.

## Solvers

Fluent is used as flow solver.
The provided mesh is triangular. When the gate bends re-meshing is performed to preserve its quality.
A script to regenerate it using Gambit is included. This script allows to change the resolution and geometrical parameters.

The structure solver is Abaqus.
The Abaqus case is built when setting up the case starting from the file *`mesh_turek.inp`* containing nodes and elements. 
This is done by running Abaqus with the *`make_inp.py`* Python script to set all parameters, such as surface definitions, material parameters, boundary conditions and time step information.
The result of the setup is a completed input file *`case_turek.inp`*.
The Abaqus element type used is CPE8R.

To exchange information between the solvers on the fluid-structure interface, the use of mappers is required.
In the structure solver wrapper, a linear interpolation mapper is used to interpolate in the x- and y-direction from and to the coupled solver.

## References
<a id="1">[1]</a> 
[Turek S., Hron J., "Proposal for Numerical Benchmarking of Fluid-Structure Interaction between an Elastic Object and Laminar Incompressible Flow", Bungartz HJ., Schäfer M. (eds) Fluid-Structure Interaction. Lecture Notes in Computational Science and Engineering, vol. 53. Springer, Berlin, Heidelberg, 2006.](https://doi.org/10.1007/3-540-34596-5_15)