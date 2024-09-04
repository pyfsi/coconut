# Yarn case with FluentALM and Abaqus3D - Steady

Example calculates the static yarn deformation under the action of an air jet using the ALM formulation in the flow solver and line loads in the structural solver, see the Figure below.

<img src=".velocity_xy_deformed_yarn.png" alt="Velocity contour and deformed yarn at the end of the simulation." width="750"/>

## Solvers

### Fluent

The setup of this example starts from a mesh file and calculates the steady-state flow field without the presence of the yarn as initial condition to the FSI simulation. 
Only upon running the actual coupled simulation, the yarn is included as actuator line.
Since the mesh counts about 1M cells, this example is hardcoded to run on a machine with at least 40 cores.

### Abaqus

The structural model is at this point not yet based on the one derived in Bral et al. [[1](#1)], but will be incorporated later on.

## References
<a id="1">[1]</a> 
[Bral A., Daelemans L. and Degroote J., "A novel technique to simulate and characterize yarn mechanical behavior based on a geometrical fiber model extracted from microcomputed tomography imaging", Textile Research Journal, vol. 93(9-10), pp. 2042-2062, 2023.](https://doi.org/10.1177/00405175221137009)
