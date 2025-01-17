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


## Parameters

All parameters used in the orginal [Fluent solver wrapper](../fluent/fluent.md) remain available.
A new subdictionary with keyword `PC` should be provided, however, containing the following keywords:

|               parameter | type  | description                                                                                                                                 |
|------------------------:|:-----:|---------------------------------------------------------------------------------------------------------------------------------------------|
|       `moving_boundary` | bool  | (optional) Default: `true`. Should normally be always `true`, as phase change problems are moving boundary problems.                        |
|         `ini_condition` | float | Scalar value as initial condition for the output thermal boundary condition: temperature (in K) or heat flux (in W/m$\cdot$K).              |
|                `latent` | float | Latent heat of the phase change material (PCM).                                                                                             |
|             `melt_temp` | float | Melting temperature of the phase change material (PCM).                                                                                     |
| `end_of_setup_commands` |  str  | (optional) Fluent journal command(s) to be executed after the setup is finished, can be used to indicate the cell height for mesh layering. |


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

Within the solid domain, only the energy equation is solved with a source term corresponding to the enthalpy method for phase change [[1](#1)].
The enthalpy method ensures a correct temperature field and derived heat fluxes in the solid phase.
Consequently, an additional energy source term is required to account for the *'excess'* cell enthalpy in the interface cells inherited from the previous time step.
This *'excess'* enthalpy represents the latent part of the cell enthalpy, which is added to the cell's sensible enthalpy at melting temperature.
If not removed by a source term (which acts as a sink during melting), the cell enthalpy of boundary cells would continue to increase due to the incoming
heat flux across multiple time steps. This would lead to an increase in the liquid fraction despite the mesh updates, which should remove all molten material each time step.

The heat flux profile returned by the liquid solver serves as input for the solid solver through the *`set_heat_flux`* UDF.
The *`update_cell_enthalpy`* UDF stores the cell enthalpy of cells adjacent to the interface, which is necessary for the energy source term.
The energy equation with the correct source term is then solved, enabling the calculation of the new interface displacement.

This is done through the *`write_displacement`* UDF, where the Stefan condition is applied to each cell face at the interface.
The displacement follows from the interface velocity $v_{itf}$, which can be determined by enforcing the Stefan condition:
$$
\rho \cdot L \cdot v_{itf} = -k_L \nabla T |^L + k_S \nabla T |^S
$$
with $\rho$ the density of the PCM, $L$ the latent heat, $-k_L \nabla T|^L$ the interface heat flux at the liquid side and $-k_S \nabla T|^S$ the same at the solid side.
The face displacement is written to a file and the fluent script of the solver wrapper subsequently converts this to node displacement before it is passed on as output.
A built-in mapper in CoCoNuT facilitates the [face-to-node mapping](../../mappers/mappers.md#linearconservative).
These node displacements are then returned to the liquid solver for the next coupling iteration.

If the solid is at the melting temperature, the interface temperature is assumed to be at the melting temperature as well, eliminating the need to exchange interface temperature as a variable.
However, when the solid is subcooled below the melting temperature, it may require initial sensible heating before phase change occurs.
In such cases, the problem initially becomes a conjugate heat transfer problem, necessitating the exchange of temperature from the solid domain to the liquid domain.
The transition from sensible heating to latent heating also requires special treatment.

In summary, the coupling strategy is illustrated in the figure below.

![](images/pc_coupling.png "Coupling strategy between liquid and solid domain for phase change problems")


### Files created during simulation

In these file conventions, A is the time step number and B the Fluent thread ID.

-   Fluent case and data files are saved as files of the form *`case_timestepA.cas`* and *`case_timestepA.dat`*.
-   Current node and face coordinates are passed from Fluent to CoCoNuT with files of the form *`nodes_timestepA_threadB.dat`* and *`faces_timestepA_threadB.dat`*.
-   Face displacement is passed from Fluent to CoCoNuT with files of the form *`displacement_timestepA_threadB.dat`*.
-   The new node coordinates are passed from CoCoNuT to Fluent with files of the form *`nodes_update_timestepA_threadB.dat`*.
-   Heat flux is passed from Fluent to CoCoNuT with files of the form *`heat_flux_timestepA_threadB.dat`*.
-   Temperature is passed from Fluent to CoCoNuT with files of the form *`temperature_timestepA_threadB.dat`*.
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

