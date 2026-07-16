# Two-dimensional validation: paraffin melting

## Introduction & Physics
The second validation problem, discussed in [[1](#1)], considers a two-dimensional square cavity containing a phase change material (PCM) modelled on octadecane.
The cavity has a width and height of 0.04 m. Initially, the PCM is completely solid at 291.13 K, which is 10 K below its melting temperature.
At t = 0, the left wall is heated to 311.13 K (10 K above melting), while the right wall is kept at 291.13 K.
The top and bottom walls are perfectly insulated, and a no-slip boundary condition is applied to all walls.
To compensate for the volume expansion of the PCM during phase change, a small outlet (0.0004 m wide) at zero gauge pressure is located on the left side of the top wall.

The partitioned simulation is initialised at a global liquid fraction of 0.02.
Because only a narrow liquid gap exists near the heated wall at this stage, convective currents are negligible.
The initial heat transfer is purely conductive, resulting in a straight solid-liquid interface. 
The domain is thus divided into a rectangular liquid region (0.0008 m wide) and a solid region (0.0392 m wide),
with the initial temperature field defined by the analytical solution of the two-phase Stefan problem.

## Parameters
*   **Material properties (Octadecane)**:
    *   Solid density: 867.914 kg/m³
    *   Liquid density: $979.826 - 0.674 \cdot T$ kg/m³
    *   Solid specific heat capacity: $-1029 + 9.797 \cdot T$ J/(kg K)
    *   Liquid specific heat capacity: $3247 - 8.861 \cdot T + 1.821 \cdot 10^{-2} \cdot T^2$ J/(kg K)
    *   Solid thermal conductivity: 0.334 W/(m K)
    *   Liquid thermal conductivity: $0.246 - 3.121 \cdot 10^{-4} \cdot T$ W/(m K)
    *   Liquid dynamic viscosity: $10^{-3} \cdot \exp(-5.353 + 2026.013 / T)$ kg/(m s)
    *   Melting temperature: 301.13 K
    *   Latent heat: 236980 J/kg
*   **Domain & setup parameters**:
    *   Hot wall temperature: 311.13 K
    *   Cold wall temperature: 291.13 K
    *   Domain dimensions: 0.04 m width by 0.04 m height
    *   Time step size: 0.1 s

## Solvers
*   **Liquid domain**: ANSYS Fluent 2024R2 (pressure-based solver) is used. The domain employs an unstructured mesh with two rows of rectangular boundary layer cells at the interface. Spring-based smoothing and the unified remeshing algorithm are applied to manage volume changes, both limited to a maximum cell skewness of 0.4. A second-order upwind scheme is used for convective terms, the PISO algorithm for velocity-pressure coupling, least-squares cell-based gradients, and body-force-weighted pressure interpolation. Time stepping uses a second-order implicit scheme.
*   **Solid domain**: ANSYS Fluent 2024R2 is also employed for the solid domain. It uses a structured mesh of 800 vertical and 400 horizontal subdivisions, featuring corresponding boundary layer cells on the solid side of the interface. Diffusion-based smoothing is applied. Boundary layer nodes in both domains follow the interface displacement exactly to maintain a constant thickness.

## Algorithmic Settings
*   **Coupling algorithm**: Explicit coupling is applied for the two domains.
*   **Predictor**: Linear predictor.
*   **Convergence criterion**: Not applicable for explicit coupling.

## References

<a id="1">[1]</a>
[V. Van Riet, W. Beyne, and J. Degroote, “Partitioned approach for convective melting of constrained, subcooled solids,” in XI International Conference on Coupled Problems in Science and Engineering (COUPLED PROBLEMS 2025), Villasimius, Italy, 2025.]
