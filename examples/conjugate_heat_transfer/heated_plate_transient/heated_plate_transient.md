# Transient heated flat plate

This is a transient analysis of a heated flat plate, similar to the steady-state case discussed in [Steady heated plate](conjugate_heat_transfer/heated_plate_steady/heated_plate_steady.md).
The only difference is that the inlet velocity varies with time. Please refer to the documentation of the steady-state case for more details on the problem.

The problem is initialized at the steady-state solution for an inlet velocity of $u_0 = 2.4 m/s$.
The inlet velocity will then sinusoidally increase from 2.4 m/s to 36 m/s in a time span of 1 s, and remains fixed at 36 m/s for another 1 s.
This fast transient behaviour will challenge the coupling algorithm as the Biot number along the interface will quickly increase from values predominantly below 1
to values much larger than 1.
The increasing velocity of the hot fluid flow will improve the heat transfer between both domains and is expected to increase the interface temperature.

## Solvers

Fluent is used for both the fluid and solid domain. Only the conduction energy equation is solved in the solid domain,
while the fluid solver considers the Navier-Stokes equations as well.
Coconut interacts with Fluent through the solver wrapper _`cht_fluent.py`_, which is developed for conjugate heat transfer applications.
The exchanged variables are heat flux and temperature. The order in which they are exchanged has influence on the stability of the coupling.

The grids are provided in the setup_files directory. The grids are structured with refinement near the interface.

## Coupling algorithm

The coupling technique used is the *interface quasi-Newton algorithm with an approximation for the inverse of the Jacobian from a least-squares model* (IQN-ILS).
The reuse parameter `q` is set to 3, which means that data from the last three time steps will be used to stabilize and accelerate the convergence.

## Predictor

The initial guess in every time step is done using the linear predictor.

## Convergence criterion

Three convergence criteria have been specified:

- The number of iterations in every time step is larger than 20.
- The residual norm of the temperature is a factor $10^{-2}$ lower than the initial value.
- The residual norm of the temperature is lower than $10^{-2}$ in absolute value.
 
When either criterion is satisfied the simulation stops.