# The two-phase Stefan problem

This example solves the two-phase Stefan problem and compares the numerical result with the analytical solution.
The two-phase Stefan problem is a one-dimensional problem where melting occurs in a (semi-infinite) slab.
The slab is initially solid at a temperature $T_S$ lower than the phase change temperature $T_m$.
The solid phase heats up and starts melting by imposing a constant boundary temperature $T_L > T_m$ at $x = 0$.
The problem is called "two-phase" since the solid phase (and not just the liquid phase as in the "one-phase" problem) will also experience a temperature gradient.
This means that conduction occurs in both the solid and liquid phases.

The problem has been accomodated for numerical purposes by using a finite domain, where the wall at the opposite side is at the initial solid temperature $T_S$.
The challenge of the problem lies in the tracking of the moving interface.
A script _`validate_two_phase_Stefan.py`_ is provided to compare the results with the analytical solution [[1](#1)] and the
results from simulations in Fluent using the enthalpy method for phase change.

The PCM (phase change material) parameters are:

- conduction coefficient of liquid phase: 1.5 W/m$\cdot$K
- conduction coefficient of solid phase: 0.024 W/m$\cdot$K
- heat capacity: 2500 J/kg$\cdot$K
- density: 870 kg/m³
- latent heat of phase change: 179000 J/kg

The domain parameters are:

- domain length = 0.1 m
- hot boundary temperature $T_L$: 309.15 K 
- cold boundary temperature $T_S$: 289.15 K
- phase change temperature $T_m$: 299.15 K

The total simulated time is 7200 s (2 h) in time steps of 0.1 s.
Due to the slow physics of the problem, the example requires a lot of time steps and takes a very long time to run.

The problem is initialized with the temperature profile of the analytical solution at a liquid fraction of 50 % and wil run
to a liquid fraction of approximately 60 % in a simulation time of 2 h.

The figure below shows the liquid fraction as function of time.
The coconut simulation results are in closer agreement with the analytical solution than the results from the enthalpy method in Fluent.

![PC1](images/stefan_two_phase_liq_frac.png "Liquid fraction as a function of time for the two-phase Stefan problem.")

## Solvers

Fluent is used for both the liquid and solid phase.
Coconut interacts with Fluent through the solver wrapper _`pc_fluent.py`_, which is being developed for phase change applications.
This solver wrapper manages both Fluent solver processes for the solid and liquid domain.
In the case where the solid phase can be subcooled below phase change temperature, interface temperature should be exchanged as well besides heat flux and displacement.
The interface can no longer be assumed at phase change temperature in the initial stages of the simulation.

The grids are provided in the setup_files directory. They represent respectively 50 % (liquid) and 50 % (solid) of the domain.
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
-   The residual norm of the displacement and temperature combined is a factor $10^{-4}$ lower than the initial value.
 
When either criterion is satisfied the simulation stops.
 
## References
<a id="1">[1]</a>
[Alexiades V., Solomon, A. D., "2.2. The two-phase problem on a semi-infinite slab", in "Mathematical modeling of melting and freezing processes", London, Taylor & Francis, 1993, pp. 46-58]