# Tube case with Fluent2D and Abaqus2D - Steady

This example calculates the flow inside and the deformation and stresses of a straight flexible tube, when a steady pressure difference is applied over the tube.
This done by using Fluent and Abaqus, both on an axisymmetric case.

The test example is similar in setup to 'tube_fluent2d_abaqus2d'. The only difference is that a steady problem is simulated.
Here the most important differences and peculiarities of a steady case are highlighted.

## General settings
Although the caculation is steady, a `delta_t` is still required. Its value is abritrary and usualy 1.0 is used.
`timestep_start` is required as well and is normally equal to 0.

## Coupling algorithm

The coupling technique used is the *interface quasi-Newton algorithm with an approximation for the inverse of the Jacobian from a least-squares model* (IQNI-LS).
The paramter `q` is not used no as there is only one time step. Note that in a steady caculation, the models `ls` (IQN-ILS) and `mv` (IQN-MVJ) are identical.

## Predictors

A predictor is still required, but not used as only one time step is calculated.

## Solvers

Of course the suplied case files in both Fluent and Abaqus also need to be steady. 
In Abaqus this can be done using a *Static* step, typically with subcycling and and a linearly ramped load.
The parameters `subcycling` and `ramp` are set on true in the json-file.
As such Abaqus performs subiterations in each coupling iterations in which the load is increased linearly over the step.
The parameters `minInc`, `initialInc`, `maxNumInc` and `maxInc` are used to determine its behaviour.
The ramping does not occur in Abaqus itself as, amplitude references are ignored for nonuniform loads given by user subroutine DLOAD in an Abaqus/Standard analysis.
Instead, the ramping is implemented in the DLOAD subroutine itself.
For the first iteration of the first time step an initial load is required which is set to zero in the Abaqus wrapper.