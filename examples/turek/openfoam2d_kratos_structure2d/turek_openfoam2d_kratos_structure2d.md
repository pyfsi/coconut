# Turek benchmark with OpenFOAM2D and KratosStructure2D

This example calculates the well-known Turek benchmark, which consist of laminar incompressible flow around a flexible flag attached to a rigid cylinder in a 2D channel.
The material parameters and boundary conditions correspond to the FSI2 setup, as detailed by Turek and Hron [[1](#1)].
The challenging aspect of this benchmark lies in the large self-induced oscillation of the flag.
The used solvers are OpenFOAM and Kratos Multiphysics Structural Mechanics application.
Furthermore, the script _`animate_beam.py`_ can be used to visualize the interface data throughout the calculation.
A video made with ParaView is shown below.

<iframe width="560" height="315" src="https://www.youtube.com/embed/YAe1HrT89k4" title="FSI simulation of Hron &amp; Turek benchmark" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" allowfullscreen></iframe>

With these solver, the run time is significantly smaller than that of the same test case with [Fluent and Abaqus](../fluent2d_abaqus2d/turek_fluent2d_abaqus2d.md).
For details on boundary conditions and parameters refer to the example with [Fluent and Abaqus](../fluent2d_abaqus2d/turek_fluent2d_abaqus2d.md).

The accuracy of the simulation framework can be assessed by comparing with the results in Turek and Hron [[1](#1)].
A script _`evaluate_benchmark.py`_ is provided to extract the necessary quantities, which are given below using the format "mean value +/- amplitude [frequency]".

|                                  | x displacement                               | y displacement                              |
|---------------------------------:|:---------------------------------------------|---------------------------------------------|
| CoCoNuT with OpenFOAM and Kratos | <nobr>-1.4627e-02 +/-1.2248e-02 [3.9]</nobr> | <nobr>1.2494e-03 +/-8.0609e-02 [1.9]</nobr> |
|                        Benchmark | -1.4580e-02 +/-1.2440e-02 [3.8]              | 1.2300e-03 +/-8.0600e-02 [2.0]              |

|                                  | drag                                        | lift                                         |
|---------------------------------:|:--------------------------------------------|----------------------------------------------|
| CoCoNuT with OpenFOAM and Kratos | <nobr>2.1493e+02 +/-7.4956e+01 [3.9]</nobr> | <nobr>-6.7446e-01 +/-2.3692e+02 [1.9]</nobr> |
|                        Benchmark | 2.0883e+02 +/-7.3750e+01 [3.8]              | 8.8000e-01 +/-2.3420e+02 [2.0]               |


## Coupling algorithm

The coupling technique used is the *interface quasi-Newton algorithm with an approximation for the inverse of the Jacobian from a least-squares model* (IQN-ILSM).
The reuse parameter `q` is set to 10, which means that data from the last ten time steps will be used to stabilize and accelerate the convergence.

## Predictor

The initial guess in every time step is done using the quadratic predictor.

## Convergence criterion

Two convergence criteria have been specified:

-   The number of iterations in every time step is larger than 15.
-   The absolute norm of the displacement is lower than $10^{-7}$.
 
When either criterion is satisfied the simulation stops.

## Solvers

OpenFOAM is used as flow solver.
The provided mesh is quadrilateral.
Some finer mesh settings are provided in the blockMeshDict.

The structural solver is Kratos Multiphysics and uses first-order, quadrilateral, plane strain elements.

To exchange information between the solvers on the fluid-structure interface, the use of mappers is required.
In the structural solver wrapper, a linear interpolation mapper is used to interpolate in the x- and y-direction from and to the coupled solver.

## References
<a id="1">[1]</a> 
[Turek S., Hron J., "Proposal for Numerical Benchmarking of Fluid-Structure Interaction between an Elastic Object and Laminar Incompressible Flow", Bungartz HJ., Schäfer M. (eds) Fluid-Structure Interaction. Lecture Notes in Computational Science and Engineering, vol. 53. Springer, Berlin, Heidelberg, 2006.](https://doi.org/10.1007/3-540-34596-5_15)