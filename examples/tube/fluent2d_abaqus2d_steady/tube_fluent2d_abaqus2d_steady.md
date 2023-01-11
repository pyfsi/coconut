# Tube case with Fluent2D and Abaqus2D - Steady

This example calculates the flow inside and the deformation and stresses of a straight flexible tube, when a steady pressure difference is applied over the tube.
This done by using Fluent and Abaqus, both on an axisymmetric case.

The test example is similar in setup to *`tube_fluent2d_abaqus2d`*. The only difference is that a steady problem is simulated.
Here the most important differences and peculiarities of a steady case are highlighted.

## General settings
Although the calculation is steady, a `delta_t` is still required. Its value is arbitrary and usually 1.0 is used.
`timestep_start` is required as well and is normally equal to 0.

## Coupling algorithm

The coupling technique used is the *interface quasi-Newton algorithm with an approximation for the inverse of the Jacobian from a least-squares model* (IQNI-LS).
The parameter `q` is not used no as there is only one time step. Note that in a steady calculation, the models `ls` (IQN-ILS) and `mv` (IQN-MVJ) are identical.

## Predictor

A predictor is still required, but not used as only one time step is calculated.

## Solvers

Of course the supplied case files in both Fluent and Abaqus also need to be steady. 
In Abaqus this can be done using a *Static* step, typically with automatic time incrementation (subcycling) and a linearly ramped load.
The parameter `ramp` is set to `true` in the json-file. As such Abaqus performs subiterations in each coupling iterations in which the load is increased linearly over the step.
The ramping does not occur in Abaqus itself as, amplitude references are ignored for nonuniform loads given by user subroutine DLOAD in an Abaqus/Standard analysis.
Instead, the ramping is implemented in the DLOAD subroutine itself.
For the first iteration of the first time step an initial load is required which is set to zero in the Abaqus wrapper.