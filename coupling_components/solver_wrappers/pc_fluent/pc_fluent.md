# Phase change in Fluent

This is the documentation for the Fluent solver wrapper adapted for phase change simulations within CoCoNuT.
The functioning of this solver wrapper is entirely based on the original [Fluent solver wrapper](../fluent/fluent.md)
and builds further on the (thermal) features 

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
|                `latent` | float | Latent heat of the phase change material.                                                                                                   |
|             `melt_temp` | float | Melting temperature of the phase change material.                                                                                           |
| `end_of_setup_commands` |  str  | (optional) Fluent journal command(s) to be executed after the setup is finished, can be used to indicate the cell height for mesh layering. |


## Overview of operation

The solver wrapper consists of 3 files (with X the Fluent version, e.g. "v2023R1"):

-   *`X.py`*: defines the `SolverWrapperPCFluentX` class, 
-   *`X.jou`*: Fluent journal file to interactively run the phase change simulation, written in Scheme, 
-   *`udf_thermal.c`*: Fluent UDF file that implements the source terms and additional functionality (I/O) used in Fluent, written in C.


### Added functionality

The primary goal of the added functionality is to be able to exchange heat flux and interface displacement between the solvers.
In this, the heat flux is calculated at the liquid side and enters the solid side during melting.
As such, the incoming heat flux causes the solid domain to shrink and the interface to move inwards towards the solid domain.
The two largest challenges are how to calculate the interface displacement and how to accommodate the volume change in the two separate domains.
Fluent udf's have been developed to handle these challenges for both the liquid and solid domain.

The liquid side receives the interface displacement as input during phase change. The mesh is deformed accordingly with the existing *`move_nodes`* udf.
The *`set_adjacent`* udf then flags the cells at the interface and stores the new coordinates of the face nodes at the coupling interface.
These new node positions are compared to the node positions of the previous mesh update in the *`calc_volume_change`* udf to determine the swept volume of each face at the coupling interface.
The swept volume is apointed the according cells at the interface and are used for the determination of the mass and energy source terms in the flagged cells at the interface.
The flow and energy equations with the correct source terms are subsequently solved in the liquid domain and the resulting heat flux at the interface is written to a file by the *`store_heat_flux`* udf.

In the solid domai, only the energy equation with a source term for the enthalpy method is solved.
Using the enthalpy method guarantees a correct temperature field (and derived heat fluxes) in the solid phase.
Consequently, the solid solver requires an additional energy source term to remove the *'excess'* cell enthalpy in the cells adjacent to the interface from the previous time step.
This *'excess'* enthalpy is the latent part of the cell enthalpy which adds to the sensible enthalpy at melting temperature.
If not removed by a source term (which acts as a sink during melting), the cell enthalpy of the boundary cells keeps increasing due to the incoming heat flux.
This would cause the liquid fraction in the solid domain to increase despite the mesh updates, which should remove all molten material.

The heat flux profile returned by the liquid solver is used as input for the solid solver through the *`set_heat_flux`* udf.
The cell enthalpy of the cells adjacent to the interface, necessary for the energy source term, are stored by the *`update_cell_enthalpy`* udf.
The energy equation with correct source term is once again solved, which allows the new interface displacement to be calculated.

This is done through the *`write_displacement`* udf where the Stefan condition is applied to each cell face at the interface.
The face displacement, written to a file, are converted to node displacement in the fluent script of the solver wrapper.
An instance of the built-in mappers in CoCoNuT is used for this. These node displacements are then returned to the liquid solver for the next iteration.

![](images/pc_coupling.png "Coupling strategy between liquid and solid domain for phase change problems")

## Version specific documentation

### v2023R1 (23.1.0)

Base version.

### v2024R1 (24.1.0)

No changes.


