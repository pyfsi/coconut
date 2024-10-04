# The one-phase Stefan problem

This example solves the one-phase Stefan problem and compares the numerical result with the analytical solution.
The one-phase Stefan problem is a one-dimensional problem where melting occurs in a (semi-infinite) slab.
The slab is initially solid at the phase change temperature $T_m$ and starts melting by imposing a constant boundary temperature $T_L > T_m$ at $x = 0$.
The problem is called "one-phase" because only the melt experiences a temperature gradient, and as a result conduction heat transfer,
while the solid phase is at constant temperature.

The problem has been accomodated for numerical purposes by using a finite domain, where the wall at the opposite side is at the phase change temperature $T_m$.
The challenge of the problem lies in the tracking of the moving interface.
A script _`validate_one_phase_Stefan.py`_ is provided to compare the results with the analytical solution and the results from simultions in Fluent using the enthalpy-method.