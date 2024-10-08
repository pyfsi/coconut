# Phase change in Fluent

This is the documentation for the Fluent solver wrapper adapted for phase change simulations within CoCoNuT.


## How solid-liquid phase change is simulated with a partitioned approach

During phase change, two phases are present: a solid- and liquid phase.
In a partitioned approach, the phases are solved separately and interact each other at the phase change front.
For example, during melting of ice, heat flux enters from the water to the ice, which causes mass transfer of molten ice from the solid domain to the liquid domain.
Within the framework of CoCoNuT, Fluent is used as solver for both the solid and liquid domain.
In the solid domain, only the heat equation is solved while the flow equations are solved as well in the liquid domain.
The solver wrapper *pc_fluent*, however, functions differently depending on the phase.
In case of the solid domain, an incoming heat flux leads to a displacement of the interface while the liquid domain receives the interface displacement,
updates the mesh, and recalculates the outgoing heat flux.
In essence, heat flux and displacement are passed on as variables between the solvers in the manner shown in the figure below.