# Phase change in Fluent

This is the documentation for the Fluent solver wrapper adapted for phase change simulations within CoCoNuT.
The functioning of this solver wrapper is entirely based on the original [Fluent solver wrapper](../fluent/fluent.md)
and builds further on the (thermal) features introduced in the [Conjugate heat transfer solver wrapper for Fluent](../cht_fluent/cht_fluent.md).

## How solid-liquid phase change is simulated with a partitioned approach

During phase transitions, two distinct phases coexist: solid and liquid.
In a partitioned approach, these phases are modeled independently, with interactions occurring at the phase change interface.
For instance, during ice melting, heat flux from the liquid phase (water) induces mass transfer from the solid (ice) to the liquid domain.

Within the CoCoNuT framework, Fluent is employed as the solver for both solid and liquid domains.
While the heat equation is solely solved in the solid domain, both heat and flow equations are considered in the liquid domain.
However, the *pc_fluent* solver wrapper operates differently depending on the phase.

In the solid domain, an incoming heat flux results in interface displacement.
Conversely, the liquid domain receives this displacement, updates the mesh accordingly, and recalculates the outgoing heat flux.
Essentially, heat flux and interface displacement are exchanged as variables between the solvers.


## Known limitations and untested conditions

* Only 2D cases. Axisymmetric or 3D cases are currently not supported.
* Currently, only melting is supported, no solidification.
* Currently, only constrained melting cases are possible. The goal is to add mechanical coupling between the domains in the close future to allow unconstrained melting cases.
* Subcooled solids are possible, but should be initialised in such a way that the phase change interface is already at melting temperature. Heating of the solid without melting is a conjugate heat transfer problem and not possible yet.


## Parameters

All parameters used in the orginal [Fluent solver wrapper](../fluent/fluent.md) remain available.
A new subdictionary with keyword `PC` should be provided, however, containing the following keywords:

|               parameter | type  | description                                                                                                                                        |
|------------------------:|:-----:|----------------------------------------------------------------------------------------------------------------------------------------------------|
|       `moving_boundary` | bool  | (optional) Default: `true`. Should normally be always `true`, as phase change problems are moving boundary problems.                               |
|         `ini_condition` | float | Scalar value as initial condition for the output thermal boundary condition: temperature (in K) or heat flux (in W/m$\cdot$K).                     |
|                `latent` | float | Latent heat of the phase change material (PCM).                                                                                                    |
|             `melt_temp` | float | Melting temperature of the phase change material (PCM).                                                                                            |
|         `melt_enthalpy` | float | (optional) Default: enthalpy is calculated with the assumption of a constant $c_p$ value for each phase. Sensible enthalpy at melting temperature. |
|         `solid_density` | float | (optional) Default: solid density is assumed equal to liquid density. Density of the solid phase at melting temperature.                           |
| `end_of_setup_commands` |  str  | (optional) Fluent journal command(s) to be executed after the setup is finished, can be used to indicate the cell height for mesh layering.        |
|           `f2n_mapping` | dict  | Settings for the [face-to-node mapper](../../mappers/mappers.md#linearconservative).                                                               |


## Overview of operation

The solver wrapper consists of 3 files (with X the Fluent version, e.g. "v2023R1"):

-   *`X.py`*: defines the `SolverWrapperPCFluentX` class, 
-   *`X.jou`*: Fluent journal file to interactively run the phase change simulation, written in Scheme, 
-   *`udf_thermal.c`*: Fluent UDF file that implements the source terms and additional functionality (I/O) used in Fluent, written in C.


### Added functionality

The primary objective of the added functionality is to facilitate the exchange of heat flux and interface displacement between the solvers.
In this context, heat flux is calculated at the liquid side of the interface, inducing a shrinkage of the solid domain and an inward movement of the interface during melting.

The two primary challenges lie in accurately calculating interface displacement and accommodating volume changes in the two separate domains.
Fluent UDFs have been developed to address these issues within both the liquid and solid domains.

During phase change, the liquid domain receives the interface displacement as input coming from the solid solver.
The mesh is deformed accordingly using the existing *`move_nodes`* UDF.
The *`set_adjacent`* UDF then flags the cells at the interface and stores the new coordinates of face nodes at the coupling interface.
These new node positions are compared to those of the previous mesh update in the *`calc_volume_change`* UDF to determine the swept volume of each face  at the coupling interface.
The swept volume is assigned to the corresponding interface cells and used to calculate the mass and energy source terms within these cells.
The flow and energy equations, incorporating these source terms, are then solved in the liquid domain, and the resulting interface heat flux is recorded by the *`store_heat_flux`* UDF.

Within the solid domain, the conservation equations for mass, momentum and energy are solved.
The mass- and momentum equation, however, have a trivial solution, but should be included to account for the ALE terms due to mesh motion.
A temperature boundary condition is set to the malting temperature at the phase change interface, just as in the liquid domain.
The mesh is deformed and the volume change is accounted for in the same way as in the liquid domain:
using mass and source terms based on the swept volume of each cell face at the interface.

The heat flux profile returned by the liquid solver serves as input for the solid solver through the *`read_liquid_hf`* UDF.
The heat flux values are stored in the UDM (user-defined memory) of the cells adjacent to the cell faces at the interface.
These will later be used to calculate the interface displacment.

This is done through the *`write_displacement`* UDF, where the Stefan condition is applied to each cell face at the interface.
The displacement follows from the interface velocity $v_{itf}$, which can be determined by enforcing the Stefan condition:
$$
\rho \cdot L \cdot v_{itf} = -k_L \nabla T |^L + k_S \nabla T |^S
$$
with $\rho$ the density of the PCM, $L$ the latent heat, $-k_L \nabla T|^L$ the interface heat flux at the liquid side and $-k_S \nabla T|^S$ the same at the solid side.
The liquid heat flux is stored in the UDM of the cell adjacent to the face and the solid heat flux follows from the local temperature gradient in the solid domain.
The face displacement is written to a file and the fluent script of the solver wrapper subsequently converts this to node displacement before it is passed on as output.
A built-in mapper in CoCoNuT facilitates the [face-to-node mapping](../../mappers/mappers.md#linearconservative).
These node displacements are then returned to the liquid solver for the next coupling iteration.

The interface is always assumed to be at melting temperature, eliminating the need to exchange the interface temperature as a variable between the solvers.
As a result, the current implementation in CoCoNuT requires the problem to be initialised with both the solid and liquid domains already present
and the interface at melting temperature.
For example, for constrained melting in a cavity heated from one side,
this can be done by creating a small initial liquid domain and a complementary solid domain as if the liquid fraction were 0.01.
The temperature field in both domains can be initialised using the analytical solution to the Stefan problem.
This is because conductive heat transfer is dominant in the early stages of melting anyway.

Consequently, the simulation cannot start from a fully solid domain, or with an interface below melting temperature.
The latter would require exchanging temperature as a variable, as this is initially a conjugate heat transfer problem.

In summary, the coupling strategy is illustrated in the figure below.

![](images/pc_coupling.png "Coupling strategy between liquid and solid domain for phase change problems")


### Files created during simulation

In these file conventions, A is the time step number and B the Fluent thread ID.

-   Fluent case and data files are saved as files of the form *`case_timestepA.cas`* and *`case_timestepA.dat`*.
-   Current node and face coordinates are passed from Fluent to CoCoNuT with files of the form *`nodes_timestepA_threadB.dat`* and *`faces_timestepA_threadB.dat`*.
-   Face displacement is passed from Fluent to CoCoNuT with files of the form *`displacement_timestepA_threadB.dat`*.
-   The new node coordinates are passed from CoCoNuT to Fluent with files of the form *`nodes_update_timestepA_threadB.dat`*.
-   Heat flux is passed from Fluent to CoCoNuT with files of the form *`heat_flux_timestepA_threadB.dat`*.
-   Files with extension *`.coco`* are used to exchange messages between CoCoNuT and Fluent. 


## Setting up a new case

Following items should be set up and saved in the Fluent case file (this list may be non-exhaustive):

-   additional UDFs must be configured, 
-   steady/unsteady (should match with the `unsteady` parameter),
-   2D, 3D or axisymmetric (should match with the `dimensions` parameter), currently only 2D is supported,
-   dynamic mesh zones for all deforming surfaces, except for the coupling interfaces,
-   dynamic mesh smoothing parameters,
-   boundary conditions, material properties, numerical models, discretization schemes, operating conditions, turbulence modeling, convergence criteria.

A data file should also be present with the fields either initialized or containing the results of a previous calculation.

Following items are taken care of by CoCoNuT, and must therefore not be included in the Fluent case file:

-   dynamic mesh zones for the FSI interfaces (these are defined in `thread_names`),
-   the time step (`delta_t`).

## Version specific documentation

### v2023R1 (23.1.0)

Base version.

### v2024R1 (24.1.0)

No changes.

### v2024R2 (24.2.0)

No changes in operation. The Scalar_Reconstruction macro required an additional argument in the write_displacement udf.

## References

<a id="1">[1]</a>
[Voller V. R. and Prakash C., "A fixed grid numerical modelling methodology for convection-diffusion mushy region phase-change problems", Int. J. Heat Mass Transfer., vol. 30(8), pp. 1709-1719, 1987.]

