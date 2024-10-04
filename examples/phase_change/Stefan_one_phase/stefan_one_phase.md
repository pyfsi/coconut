# The one-phase Stefan problem

This example solves the one-phase Stefan problem and compares the numerical result with the analytical solution.
The one-phase Stefan problem is a one-dimensional problem where melting occurs in a (semi-infinite) slab.
The slab is initially solid at the phase change temperature $T_m$ and starts melting by imposing a constant boundary temperature $T_L > T_m$ at $x = 0$.
The problem is called "one-phase" because only the melt experiences a temperature gradient, and as a result, only one phase conducts heat
while the solid phase is at constant temperature.

The problem has been accomodated for numerical purposes by using a finite domain, where the wall at the opposite side is at the phase change temperature $T_m$.
The challenge of the problem lies in the tracking of the moving interface.
A script _`validate_one_phase_Stefan.py`_ is provided to compare the results with the analytical solution [[1](#1)] and the
results from simulations in Fluent using the enthalpy method for phase change.

The PCM (phase change material) parameters are:

- conduction coefficient of liquid phase: 1.5 W/mK
- heat capacity: 2500 J/kg$\cdot$K
- density: 870 kg/m³
- latent heat of phase change: 179000 J/kg

The domain parameters are:
- domain length = 0.1 m
- boundary temperature $T_L$: 309.15 K
- phase change temperature $T_m$: 299.15 K

The total simulated time is 36000 s (10 h) in time steps of 0.1 s.
Due to the slow physics of the problem, the example requires a lot of time steps and takes a very long time to run.

The problem is initialized with the temperature profile of the analytical solution at a liquid fraction of 1 % and wil run
to a liquid fraction of approximately 81 % in a simulation time of 10 h.

The figure below shows the liquid fraction as function of time. The numerical results agrees nearly perfectly with the analytical solution.
![PC1](images/stefan_one_phase_liq_frac.png "Liquid fraction as a function of time for the one-phase Stefan problem.")

## Solvers

Fluent is used for both the liquid and solid phase (although the solid phase is just a constant temperature field in this example).
Coconut interacts with Fluent through the solver wrapper _`pc_fluent.py`_, which is being developed for phase change applications.
This solver wrapper manages both Fluent solver processes for the solid and liquid domain.
In the special case of the solid being at phase change temperature, only heat flux and displacement are exchanged across the interface.
The interface is assumed at phase change temperature.

The grids are provided in the setup_files directory. They represent respectively 1 % (liquid) and 99 % (solid) of the domain.
As a result, thay are rectangularly shaped and contain quadrilateral cells in a structured grid.
The "end_of_setup_command" present in the _`parameters.json`_ file is necessary for the layering approach used for the dynamic mesh.

## Coupling algorithm

The coupling technique used is the *Aitken relaxation method*.
The maximal relaxation factor `omega_max` is set to 0.8.

## Predictor

The initial guess in every time step is done using the linear predictor.

## Convergence criterion

Two convergence criteria have been specified:

-   The number of iterations in every time step is larger than 10.
-   The residual norm of the displacement is a factor $10^{-4}$ lower than the initial value.
 
When either criterion is satisfied the simulation stops.
 
## References
<a id="1">[1]</a>
[Alexiades V., Solomon, A. D., "2.1. The one-phase Stefan problem", in "Mathematical modeling of melting and freezing processes", London, Taylor & Francis, 1993, pp. 33-46]