# Phase change in Fluent: Liquid Solver

This is the documentation for the Fluent liquid solver wrapper adapted for phase change simulations within CoCoNuT.
The functioning of this solver wrapper is entirely based on the original [Fluent solver wrapper](../fluent/fluent.md)
and builds further on the thermal features introduced in the [Conjugate heat transfer solver wrapper for Fluent](../cht_fluent/cht_fluent.md).
It is designed to work in tandem with the [Fluent solid solver wrapper](../pc_fluent_solid/pc_fluent_solid.md).

## How solid-liquid phase change is simulated

During phase transitions, two distinct phases coexist: solid and liquid.
In a partitioned approach, these phases are modelled independently with interactions occurring at the phase change interface.
During melting, heat flux from the liquid phase induces mass transfer from the solid to the liquid domain.

Within the CoCoNuT framework, Fluent is employed as the solver for both domains.
The *pc_fluent_liquid* solver wrapper is responsible for solving both the flow and energy equations in the fluid region.
The liquid solver receives the interface displacement calculated by the solid solver.
The mesh is deformed accordingly and the liquid solver computes the temperature and flow field to recalculate the outgoing heat flux.
Essentially, heat flux and interface displacement are exchanged as variables between the solvers.

The interface is always assumed to be at melting temperature.
This eliminates the need to exchange the interface temperature as a variable between the solvers.
As a result, the current implementation in CoCoNuT requires the problem to be initialised with both the solid and liquid domains already present and the interface exactly at melting temperature.
For example, for constrained melting in a cavity heated from one side, this can be done by creating a small initial liquid domain and a complementary solid domain as if the liquid fraction were 0.01.
The temperature field in both domains can be initialised using the analytical solution to the Stefan problem.
Consequently, the simulation cannot start from a fully solid domain or with an interface below melting temperature.

The detailed methodology of this partitioned approach for both saturated and subcooled solids is explained in [[1](#1)] and [[2](#2)], respectively.
This solver wrapper serves as the practical implementation of that theory and was used to conduct the simulations presented in those publications.

![](images/pc_coupling.png "Coupling strategy between liquid and solid domain for phase change problems")

## Known limitations and untested conditions

* Only 2D cases. Axisymmetric or 3D cases are currently not supported.
* Currently, only melting is supported. Solidification is not supported.
* Overset meshing is not supported for constrained melting. For unconstrained melting, however, it is required in the [`pc_fluent_liquid_rb`](../pc_fluent_liquid_rb/pc_fluent_liquid_rb.md) wrapper.
* Constrained melting only. This solver wrapper does not account for rigid body motion. The goal is to use the [`pc_fluent_liquid_rb`](../pc_fluent_liquid_rb/pc_fluent_liquid_rb.md) wrapper for unconstrained melting cases.
* Subcooled solids are possible but should be initialised in such a way that the phase change interface is already at melting temperature. Heating of the solid without melting is a conjugate heat transfer problem and not possible yet.
* The displacement calculated in the solid solver is used without modification in the next coupling iteration or time step while the liquid solver receives modified displacement values. Each coupling iteration modifies the liquid solver input according to the chosen coupled solver. As such, the interface position in the solid and liquid solver can slightly differ.

## Parameters

All parameters used in the original [Fluent solver wrapper](../fluent/fluent.md) remain available.
A new subdictionary with keyword `PC` should be provided containing the following keywords:

|                parameter | type  | description                                                                                                                                                                |
|-------------------------:|:-----:|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|          `ini_condition` | float | Scalar value as initial condition for the output thermal boundary condition: temperature (in K) or heat flux (in W/m$\cdot$K).                                             |
|              `melt_temp` | float | Melting temperature of the phase change material (PCM).                                                                                                                    |
|                 `latent` | float | Latent heat of the phase change material (PCM).                                                                                                                            |
|               `pcm_name` | str   | (optional) Default: **'pcm'**. The material name of the phase change material. This is required for multiphase simulations to correctly identify the material in Fluent.   |
|          `melt_enthalpy` | float | (optional) Default: `0.0`. Sensible enthalpy of the liquid phase at melting temperature. Required if using temperature-dependent specific heat ($c_p$).                    |
|          `solid_density` | float | (optional) Default: `0.0` (assumes equal density). Density of the solid phase at melting temperature, used to calculate volume change during phase change.                 |
|  `end_of_setup_commands` | str   | (optional) Fluent journal command(s) to be executed after the setup is finished. Currently used to indicate the cell height for mesh layering.                             |

## Overview of operation

The solver wrapper consists of 3 files (with X the Fluent version, e.g. "v2024R2"):

*   *`X.py`*: defines the `SolverWrapperPCFluentLiquid` class.
*   *`X.jou`*: Fluent journal file to interactively run the simulation, written in Scheme.
*   *`udf_thermal.c`*: Fluent UDF file that implements the mass and energy source terms, mesh motion, and additional I/O functionality.

### Functionality

The primary objective of the added functionality is to facilitate the exchange of heat flux and interface displacement between the solvers.
Fluent UDFs have been developed to accurately calculate interface displacement and accommodate volume changes.

During phase change, the liquid domain receives the interface displacement as input coming from the solid solver.
The mesh is deformed accordingly using the existing *`move_nodes`* UDF.
The *`set_adjacent`* UDF then flags the cells at the interface, while the *`calc_volume_change`* UDF then stores the new coordinates of the face nodes at the coupling interface.
Subsequently, these new node positions are compared to those of the previous mesh update to determine the swept volume of each face at the coupling interface.

The swept volume is assigned to the corresponding interface cells and used to calculate the mass and energy source terms within these cells.
The mass and energy source terms actively add or remove PCM mass and enthalpy to compensate for the volume change of the domain.
The flow and energy equations incorporating these source terms are then solved in the liquid domain.
The resulting interface heat flux is recorded by the *`store_heat_flux`* UDF and passed as output to the solid solver.

### Files created during simulation

In these file conventions, A is the time step number and B the Fluent thread ID.

*   Fluent case and data files are saved as *`case_timestepA.cas.h5`* and *`case_timestepA.dat.h5`*.
*   Node and face coordinates are passed from Fluent to CoCoNuT as *`nodes_timestepA_threadB.dat`* and *`faces_timestepA_threadB.dat`*.
*   The new node coordinates are passed from CoCoNuT to Fluent as *`nodes_update_timestepA_threadB.dat`*.
*   Heat flux is passed from Fluent to CoCoNuT as *`heat_flux_timestepA_threadB.dat`*.
*   Files with extension *`.coco`* are used to exchange messages between CoCoNuT and Fluent.

## Setting up a new case

The following items should be configured in the Fluent case file:

*   Additional UDFs must be compiled and hooked.
*   Steady/unsteady settings (must match the `unsteady` parameter).
*   2D planar solver. Axisymmetric and 3D are currently not supported.
*   Multiphase VOF model enabled if `multiphase` is set to `true`. 
*   Only laminar or kw-sst viscous models are supported.
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
