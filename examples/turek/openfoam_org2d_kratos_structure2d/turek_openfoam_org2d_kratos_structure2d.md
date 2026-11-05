# Turek benchmark with OpenFOAM2D and KratosStructure2D

This example calculates the well-known Turek benchmark, which consist of laminar incompressible flow around a flexible flag attached to a rigid cylinder in a 2D channel.
The material parameters and boundary conditions correspond to the FSI2 setup, as detailed by Turek and Hron [[1](#1)].
The challenging aspect of this benchmark lies in the large self-induced oscillation of the flag.
The used solvers are OpenFOAM and Kratos Multiphysics Structural Mechanics application.
Furthermore, the script _`animate_beam.py`_ can be used to visualize the interface data throughout the calculation.
A video made with ParaView is shown below.

<iframe width="560" height="315" src="https://www.youtube.com/embed/YAe1HrT89k4" title="FSI simulation of Hron &amp; Turek benchmark" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" allowfullscreen></iframe>

With these solvers, the run time is significantly smaller than that of the same test case with [Fluent and Abaqus](../fluent2d_abaqus2d/turek_fluent2d_abaqus2d.md).
For details on boundary conditions and parameters refer to the example with [Fluent and Abaqus](../fluent2d_abaqus2d/turek_fluent2d_abaqus2d.md).

The accuracy of the simulation framework can be assessed by comparing with the results in Turek and Hron [[1](#1)].
A script _`evaluate_benchmark.py`_ is provided to extract the necessary quantities, which are given below using the format "mean value +/- amplitude [frequency]".

|                                  | x displacement                             | y displacement                            |
|---------------------------------:|:-------------------------------------------|-------------------------------------------|
| CoCoNuT with OpenFOAM and Kratos | <nobr>-1.463e-02 +/-1.225e-02 [3.9]</nobr> | <nobr>1.250e-03 +/-8.061e-02 [1.9]</nobr> |
|                        Benchmark | -1.458e-02 +/-1.244e-02 [3.8]              | 1.230e-03 +/-8.060e-02 [2.0]              |

|                                  | drag                                      | lift                                       |
|---------------------------------:|:------------------------------------------|--------------------------------------------|
| CoCoNuT with OpenFOAM and Kratos | <nobr>2.149e+02 +/-7.491e+01 [3.9]</nobr> | <nobr>-5.939e-01 +/-2.370e+02 [1.9]</nobr> |
|                        Benchmark | 2.0883e+02 +/-7.375e+01 [3.8]             | 8.800e-01 +/-2.342e+02 [2.0]               |


## Coupling algorithm

The coupling technique used is the *interface quasi-Newton algorithm with an approximation for the inverse of the Jacobian from a least-squares model* (IQN-ILSM).
The reuse parameter `q` is set to 10, which means that data from the last ten time steps will be used to stabilize and accelerate the convergence.

## Predictor

The initial guess in every time step is made using the quadratic predictor.

## Convergence criterion

The time step is considered converged when the solver coupling convergence criterion is satisfied for both solvers.
This criterion checks whether the solver residual is below the solver tolerance at the start of the solver call (or after one iteration for OpenFOAM).

## Solvers

OpenFOAM is used as flow solver.
The provided mesh is quadrilateral.
Some finer mesh settings are provided in the blockMeshDict.

The structural solver is Kratos Multiphysics. First-order, quadrilateral, plane strain elements are used.

To exchange information between the solvers on the fluid-structure interface, the use of mappers is required.
In the structural solver wrapper, a linear interpolation mapper is used to interpolate in the x- and y-direction from and to the coupled solver.

## References
<a id="1">[1]</a> 
[Turek S., Hron J., "Proposal for Numerical Benchmarking of Fluid-Structure Interaction between an Elastic Object and Laminar Incompressible Flow", Bungartz HJ., Schäfer M. (eds) Fluid-Structure Interaction. Lecture Notes in Computational Science and Engineering, vol. 53. Springer, Berlin, Heidelberg, 2006.](https://doi.org/10.1007/3-540-34596-5_15)