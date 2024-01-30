# Lid-driven cavity case with OpenFOAM2D and KratosStructure2D

This case is an FSI-adaptation of the famous lid-driven cavity case, calculated using OpenFOAM and KratosMultiphysics StructuralMechanicsApplication.
For more information with respect to this case refer to [this example](../fluent2d_kratos_structure2d/lid_driven_cavity_fluent2d_kratos_structure2d.md).

Reference solutions are available from Mok [[1](#1)] and Valdes [[2](#2)].
The figure belows shows a comparison with the solution of the examples with Fluent and OpenFOAM.

![comparison](images/lid_driven_cavity_comparison_openfoam.png "Comparison of y-displacement of the central point of the flexible bottom with the reference solutions")

The following figures show contour plots of the pressure and velocity for this example (with Paraview).

![velocity](images/lid_driven_cavity_velocity_openfoam.gif "Animation of velocity produced with Paraview")
![pressure](images/lid_driven_cavity_pressure_openfoam.gif "Animation of pressure produced with Paraview")

## Coupling algorithm

The coupling technique used is the *interface quasi-Newton algorithm with an approximation for the inverse of the Jacobian from a least-squares model* (IQNI-LS).
No reuse is employed, so the reuse parameter `q` is set to 0.

## Predictor

The initial guess in every time step is done using the linear predictor.

## Convergence criterion

Two convergence criteria have been specified:

-   The number of iterations in every time step is larger than 10.
-   The residual norm of the displacement is smaller than $10^{-9}$.

When either criterion is satisfied the simulation stops.

## Solvers

The flow solver is the OpenFOAM solver pimpleFOAM.
A quadrilateral mesh is used with 32 cells on each side of the cavity, made with blockMesh.

The structure solver is KratosMultiphysics StructuralMechanicsApplication (abbreviated KratosStructure).
The deformable bottom is meshed with 2 layers of 32 _TotalLagrangianElement2D4N_ elements.
The movement of the left and right side is constrained.

To exchange information between the solvers on the fluid-structure interface, the use of mappers is required.
Because a two-dimensional calculation in OpenFOAM has one cell in the depth direction, the coordinates of the nodes are actually three-dimensional.
Therefore, the displacement serving as input to the flow solver is first mapped using the mapper [MapperDepth2DTo3D](../../../coupling_components/mappers/mappers.md#mapperdepth2dto3d) and subsequently mapped using a radial basis mapper.
The resulting pressure and traction forces are located on the face centers, which all lie in the same plane as where the forces are applied in the structure solver.
This means no additional transformer is required and only a radial basis mapper has to be applied.
Mapping the flow solver instead of the structure solver, means that the coupled solver will work with the interface containing the displacements stored in Kratos nodes and not the interface which has stored the displacements in twice as many OpenFOAM nodes.


## References
<a id="1">[1]</a>
[Mok D. P., “Partitionierte Lösungsansätze in der Strukturdynamik und der Fluid−Struktur−Interaktion”, PhD thesis: Institut für Baustatik, Universität Stuttgart, 2001.](https://elib.uni-stuttgart.de/handle/11682/164)

<a id="2">[2]</a>
[Valdés G. and Gerardo J., “Nonlinear Analysis of Orthotropic Membrane and Shell Structures Including Fluid-Structure Interaction”, PhD thesis: Universitat Politècnica de Catalunya, 2007.](https://www.tdx.cat/handle/10803/6866)