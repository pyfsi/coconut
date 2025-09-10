# Nylon yarn case

This example contains the setup for the case described in Section 4.1 of Bral et al. [[1](#1)], as well as the resulting pickle file.
In the folder *`post_processing`*, one can find the scripts used to generate Figure 6 of the paper. The parameter file is set up for $\epsilon = 3R$.

## Adaptations to the flow solver

Before running the case, the correct aerodynamic force coefficients should be selected.
Therefore, change lines 597, 643 and 648 of the file *`fluent_alm/alm.c`* (in the *`solver_wrappers`* directory) to:

```C
real Re, p_dyn, c_f, c_a, c_td, multiplier, mag, vol;

// [lines 598-642]

c_f = fmax(M_PI*0.27*pow(0.5*Re, -0.61), 0.02);
c_a = fabs(cos(air_values[point][SCALAR_THETA])) * c_f;

// [lines 644-647]

c_td = pow(sin(air_values[point][SCALAR_THETA]), 2.0) * (1.18 + 3.4*pow(Re, -0.89) + 0.98*pow(Re, -0.5) - 1/(1/(0.0004*Re) + Re/1100.0)) + fabs(sin(air_values[point][SCALAR_THETA])) * c_f;
```

## References
<a id="1">[1]</a> 
Bral A., Daelemans L. and Degroote J., "Modelling the fluid-structure interactions of a hairy yarn in air-jet weaving: a multiscale approach", Under review at International Journal for Numerical Methods in Engineering, 2025.