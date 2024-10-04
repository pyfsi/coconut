# Steady-state heated flat plate

The heated flat plate is a well-known example of conjugate heat transfer with an analytical solution by Luikov [[1](#1)] to benchmark the numerical solution.
A plate of certain thickness and length has a fixed temperature boundary condition imposed on its bottom side and a fluid flow runs parallel to the top of the plate.
The fluid flow has a fixed temperature at its vertical inlet and a constant inlet velocity in the steady-state case.
The figure below depicts the situation, the relevant parameters and some of their values.

![luikov_setup](images/cht_setup.svg "Geometry and boundary conditions of the flat plate solid and fluid domain [[2](#2)].")

The boundary conditions and geometry parameters are:

- thickness b: 0.01 m
- length L: 0.2 m
- bottom temperature $T_b$: 600 K
- inlet temperature: 1000 K
- inlet velocity $u_0$: 12 m/s
- static outlet pressure: 1.03$\cdot$10$^5$ Pa

The solid and fluid material properties are

- solid conduction coefficient $\lambda_s$: 0.2876 W/m$\cdot$K
- fluid conduction coefficient: 0.06808 W/m$\cdot$K
- density: 0.3525 kg/m³
- heat capacity: 1142.6 J/kg$\cdot$K
- dynamic viscosity: 3.95$\cdot$10$^{-5}$

Only a single time step is calculated as this is a steady state problem.

The challenge in this problem follows from the evolution of the Biot number $Bi$ along the interface. At the leading edge, $Bi > 1$, while at the trailing edge of the plate,
$Bi < 1$. This leads to conflicting requirements in the solver coupling to provide stability and challenges the coupling algorithm.
For more details, the reader is referred to Van Riet et al. [[2](#2)].

A script _`validate_one_phase_Stefan.py`_ is provided to plot the steady-state temperature and heat flux distribution along the interface.
The figure below shows the steady-state temperature distribution in the y-direction through the middle of the plate (x = 0.1 m)
as achieved by the simulation and two analytical approximations.

![CHT1](images/cht_steady_temp.png "Lateral temperature distribution at x = 0.1 m.")

## Solvers

Fluent is used for both the fluid and solid domain. Only the conduction energy equation is solved in the solid domain,
while the fluid solver considers the Navier-Stokes equations as well.
Coconut interacts with Fluent through the solver wrapper _`cht_fluent.py`_, which is developed for conjugate heat transfer applications.
The exchanged variables are heat flux and temperature. The order in which they are exchanged has influence on the stability of the coupling.

The grids are provided in the setup_files directory. The grids are structured with refinement near the interface.

## Coupling algorithm

The coupling technique used is the *Aitken relaxation method*.
The maximal relaxation factor `omega_max` is set to 0.65.

## Predictor

The initial guess in every time step is done using the linear predictor.

## Convergence criterion

Two convergence criteria have been specified:

-   The number of iterations in every time step is larger than 20.
-   The residual norm of the temperature is a factor $10^{-3}$ lower than the initial value.
 
When either criterion is satisfied the simulation stops.
 
## References

<a id="1">[1]</a>
[Luikov A. v., “Conjugate convective heat transfer problems.”, Int. J. Heat Mass Transfer, 17(2), pp. 257–265, 1974.]

<a id="2">[2]</a>
[Van Riet V., Beyne W., De Paepe M., Degroote J., "Convergence behaviour of partitioned methods for conjugate heat transfer problems", in ECCOMAS 2024, Proceedings, Lisbon, Portugal, 2024.]