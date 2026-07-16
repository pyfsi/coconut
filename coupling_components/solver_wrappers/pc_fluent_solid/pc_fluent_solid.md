# Phase change in Fluent: Solid Solver

This is the documentation for the Fluent solid solver wrapper adapted for phase change simulations within CoCoNuT.
The functioning of this solver wrapper focuses entirely on solving the heat equation within the solid domain and calculating the interface regression during melting.
It is designed to work in tandem with the [Fluent liquid solver wrapper](../pc_fluent_liquid/pc_fluent_liquid.md).

## How solid-liquid phase change is simulated

During phase transitions, two distinct phases coexist: solid and liquid.
In a partitioned approach, these phases are modelled independently with interactions occurring at the phase change interface.
During melting, heat flux from the liquid phase induces mass transfer from the solid to the liquid domain.

The *pc_fluent_solid* solver wrapper receives the heat flux profile from the liquid solver.
Using the Stefan condition, it calculates the interface displacement and handles the shrinkage of the solid domain.
Because the displacement is calculated per face, this wrapper utilises an internal CoCoNuT mapper to convert face displacements into node displacements before returning them as output to the liquid solver.
Unlike the liquid solver, the solid solver does not resolve fluid flow or multiphase interactions.

The interface is always assumed to be at melting temperature.
This eliminates the need to exchange the interface temperature as a variable between the solvers.
As a result, the current implementation in CoCoNuT requires the problem to be initialised with both the solid and liquid domains already present and the interface exactly at melting temperature.
Consequently, the simulation cannot start from a fully solid domain or with an interface below melting temperature.

The detailed methodology of this partitioned approach for both saturated and subcooled solids is explained in [[1](#1)] and [[2](#2)], respectively.
This solver wrapper serves as the practical implementation of that theory and was used to conduct the simulations presented in those publications.

![](images/pc_coupling.png "Coupling strategy between liquid and solid domain for phase change problems")

## Known limitations and untested conditions

* Only 2D cases. Axisymmetric or 3D cases are currently not supported.
* Currently, only melting is supported. Solidification is not supported.
* The solid wrapper does not support multiphase simulations as these are not deemed physically possible. If a multiphase Fluent case is detected, the wrapper will raise an error.
* Subcooled solids are possible but should be initialised in such a way that the phase change interface is already at melting temperature. Heating of the solid without melting is a conjugate heat transfer problem and not possible yet.

## Parameters

All parameters used in the original [Fluent solver wrapper](../fluent/fluent.md) remain available.
A new subdictionary with keyword `PC` should be provided containing the following keywords:

|                 parameter | type  | description                                                                                                                                                    |
|--------------------------:|:-----:|----------------------------------------------------------------------------------------------------------------------------------------------------------------|
|            `ini_condition`| float | Scalar value as initial condition for the output variable (typically initial displacement, usually `0.0`).                                                     |
|                   `latent`| float | Latent heat of the phase change material (PCM). Mandatory for the solid solver to calculate the Stefan condition.                                              |
|              `melt_temp`  | float | Melting temperature of the phase change material (PCM).                                                                                                        |
|          `melt_enthalpy`  | float | (optional) Default: `0.0`. Sensible enthalpy of the solid phase at melting temperature. Required if using temperature-dependent specific heat ($c_p$).         |
|            `f2n_mapping`  | dict  | Settings for the internal [face-to-node mapper](../../mappers/mappers.md#linearconservative), required to convert calculated face regression into node motion. |
| `end_of_setup_commands`   | str   | (optional) Fluent journal command(s) to be executed after the setup is finished. Can be used to indicate the cell height for mesh layering.                    |

## Overview of operation

The solver wrapper consists of 3 files (with X the Fluent version, e.g. "v2024R2"):

*   *`X.py`*: defines the `SolverWrapperPCFluentSolid` class.
*   *`X.jou`*: Fluent journal file to interactively run the simulation, written in Scheme.
*   *`udf_thermal.c`*: Fluent UDF file that implements the Stefan condition, temperature gradients, and mesh motion.

### Added functionality

Within the solid domain, the conservation equations for mass, momentum and energy are solved.
The mass and momentum equations yield a trivial solution but are included to account for the Arbitrary Lagrangian-Eulerian (ALE) terms resulting from mesh motion.
A temperature boundary condition is set to the melting temperature at the phase change interface.
The mesh is deformed and the volume change is accounted for using mass and energy source terms based on the swept volume of each cell face at the interface.

The heat flux profile returned by the liquid solver serves as input for the solid solver through the *`read_liquid_hf`* UDF.
The heat flux values are stored in the user-defined memory (UDM) of the cells adjacent to the cell faces at the interface.
These are later used to calculate the interface displacement.

This is done through the *`write_displacement`* UDF where the Stefan condition is applied to each cell face at the interface.
The displacement follows from the interface velocity $v_{itf}$ which can be determined by enforcing the Stefan condition:
$$
\rho \cdot L \cdot v_{itf} = -k_L \nabla T |^L + k_S \nabla T |^S
$$
with $\rho$ the density of the PCM, $L$ the latent heat, $-k_L \nabla T|^L$ the interface heat flux at the liquid side and $-k_S \nabla T|^S$ the same at the solid side.
The liquid heat flux is retrieved from the UDM of the cell adjacent to the face and the solid heat flux follows from the local temperature gradient in the solid domain.

The face displacement is written to a file.
The Python wrapper reads this face-based data and utilises a built-in CoCoNuT mapper ([face-to-node mapping](../../mappers/mappers.md#linearconservative)) to translate it into a continuous node displacement field.
These mapped node displacements are then returned to the liquid solver for the next coupling iteration.

### Files created during simulation

In these file conventions, A is the time step number and B the Fluent thread ID.

*   Fluent case and data files are saved as *`case_timestepA.cas.h5`* and *`case_timestepA.dat.h5`*.
*   Node and face coordinates are passed from Fluent to CoCoNuT as *`nodes_timestepA_threadB.dat`* and *`faces_timestepA_threadB.dat`*.
*   Face displacement is extracted from Fluent via *`displacement_timestepA_threadB.dat`*.
*   The internally mapped node displacement is written to *`nodes_update_timestepA_threadB.dat`* to deform the solid mesh.
*   Files with extension *`.coco`* are used to exchange messages between CoCoNuT and Fluent.
*   Node displacement is saved for restart of the solid solver with files of the form *`displacement_restart_timestepA_threadB.dat`*.

## Setting up a new case

The following items should be configured in the Fluent case file:

*   Additional UDFs must be compiled and hooked.
*   Steady/unsteady settings (must match the `unsteady` parameter).
*   2D planar solver. Axisymmetric and 3D are currently not supported.
*   Single-phase setup. The VOF model must be disabled.
*   Dynamic mesh zones for all deforming surfaces, except for the coupling interfaces.
*   Boundary conditions, material properties, numerical models, and operating conditions.

The following items are handled by CoCoNuT and must not be included in the saved Fluent case file:

*   Dynamic mesh zones for the FSI interfaces (these are defined in `thread_names`).
*   The physical time step size (`delta_t`).

## Version specific documentation

*   **v2023R1 (23.1.0)**: Out of service.
*   **v2024R1 (24.1.0)**: Out of service.
*   **v2024R2 (24.2.0)**: No changes in operation. The Scalar_Reconstruction macro required an additional argument in the write_displacement udf.

## References

<a id="1">[1]</a>
[V. Van Riet, W. Beyne, and J. Degroote, “A partitioned interface-tracking method for convective melting of constrained saturated solids,” Int. J. Heat Mass Transfer., vol. 254, pp. 127659, 2026.]

<a id="2">[2]</a>
[V. Van Riet, W. Beyne, and J. Degroote, “Partitioned approach for convective melting of constrained, subcooled solids,” in XI International Conference on Coupled Problems in Science and Engineering (COUPLED PROBLEMS 2025), Villasimius, Italy, 2025.]
