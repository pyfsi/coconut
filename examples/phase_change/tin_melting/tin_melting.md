# Two-dimensional validation: tin melting

## Introduction & Physics
The problem of tin melting in a closed cavity is a widely used benchmark for validating convection-diffusion phase change simulations, originally proposed by Hannoun et al. [[1](#1)] and used to validate the partitioned approach in [[2](#2)].
A square cavity contains saturated tin that is initially solid at its melting temperature of 505 K.
A hot temperature of 508 K is applied to the left wall while the right wall maintains the solid temperature of 505 K.
The top and bottom walls are considered adiabatic. The Boussinesq approximation is applied to model natural convection in the liquid domain.
The melting process initially occurs as a pure one-dimensional conduction problem before multi-cellular convection starts. 
As a result, the initial condition assumes a straight interface forming a thin liquid layer of 0.001 m width.

## Parameters
*   **Material properties (Tin)**:
    *   Density: 7500 kg/m³
    *   Specific heat capacity: 200 J/(kg K)
    *   Latent heat: 60000 J/kg
    *   Thermal conductivity: 60 W/(m K)
    *   Melting temperature: 505 K
    *   Dynamic viscosity: 0.006 Ns/m²
    *   Thermal expansion coefficient: 2.67e-4 K⁻¹
*   **Domain & setup parameters**:
    *   Hot wall temperature: 508 K
    *   Domain dimensions: 0.1 m width by 0.1 m height
    *   Gravitational acceleration: 10 m/s²
    *   Initial liquid width: 0.001 m
    *   Initial global liquid fraction: 0.01
    *   Time step size: 0.01 s

## Solvers
*   **Liquid domain**: The finite volume solver ANSYS Fluent 2024R2 is employed for the liquid phase. The domain is discretised using an unstructured mesh featuring a specific boundary layer grid at the interface. This boundary layer is 0.0001 m thick and contains two rectangular cell layers to enhance data exchange and the application of source terms. To compensate for the volume increase during melting, a spring-based smoothing and unified remeshing method is applied with a maximum cell skewness of 0.7. The pressure-based solver uses a second-order upwind scheme for convective terms, the PISO algorithm for velocity and pressure coupling, a least-squares cell-based method for gradients, and PRESTO! for pressure interpolation. Time stepping is handled by a second-order implicit scheme.
*   **Solid domain**: The solid domain is fully saturated and possesses no thermal gradients. A Python-based solver is therefore used to enforce the Stefan condition at the interface.

## Algorithmic Settings
*   **Coupling algorithm**: An implicit coupling scheme using Aitken relaxation is applied.
*   **Predictor**: A linear predictor is utilised for the interface displacement.
*   **Convergence criterion**: Convergence is evaluated using a combined 'OR' condition. The coupling loop converges if the maximum number of iterations (10) is reached, or if the interface displacement meets either a relative L2 norm tolerance of 1$^{-3}$ or an absolute L2 norm tolerance of 1$^{-15}$.

## References

<a id="1">[1]</a>
[N. Hannoun, V. Alexiades, T. Z. Mai, "Reference solution for the two-dimensional melting problem," Int. J. Heat Mass Transfer, 2005.]

<a id="2">[2]</a>
[V. Van Riet, W. Beyne, and J. Degroote, “A partitioned interface-tracking method for convective melting of constrained saturated solids,” Int. J. Heat Mass Transfer., vol. 254, pp. 127659, 2026.]
